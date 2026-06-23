```python
#!/usr/bin/env python3
"""
Compare translated protein FASTA files against their original protein FASTA files.

For each translated FASTA, this script looks for the matching original FASTA,
matches records by species key, globally aligns translated/original sequences,
and writes a failure-only CSV when problems are detected.

Failures include:
    1) Missing original FASTA file.
    2) Species present in translated FASTA but missing from original FASTA.
    3) Pairwise identity below the specified threshold, or no comparable residues.

Output CSVs contain only failing rows.
"""

import argparse
import csv
import os
import re
from typing import Dict, Tuple, List

from Bio import SeqIO
from Bio.Align import PairwiseAligner


def clean_seq(s: str) -> str:
    """
    Normalize a sequence string for comparison.

    Removes whitespace, converts to uppercase, and drops a terminal stop codon
    if the sequence ends with '*'.
    """
    s = re.sub(r"\s+", "", s).upper()
    # drop terminal stop
    if s.endswith("*"):
        s = s[:-1]
    return s


def species_key_before_first_underscore(rec) -> str:
    """
    Extract the species key from a FASTA record ID.

    Uses the first whitespace-delimited token of the record ID, then keeps
    everything before the first underscore.
    """
    token = rec.id.split()[0]
    return token.split("_", 1)[0]


def load_by_species_key(path: str) -> Dict[str, Tuple[str, str]]:
    """
    Load a FASTA file into a dictionary keyed by species key.

    Returns:
        Dict mapping species_key -> (record_id, cleaned_sequence)

    If multiple records have the same species key, only the first one
    encountered is kept.
    """
    #### returns species_key -> (record_id, sequence)
    #### keeps the first record encountered for a species_key
    d: Dict[str, Tuple[str, str]] = {}
    for rec in SeqIO.parse(path, "fasta"):
        key = species_key_before_first_underscore(rec)
        if key not in d:
            d[key] = (rec.id, clean_seq(str(rec.seq)))
    return d


def build_aligner() -> PairwiseAligner:
    """
    Build and configure the global pairwise aligner used for sequence comparison.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    return aligner


def align_to_gapped_strings(aligner: PairwiseAligner, s1: str, s2: str) -> Tuple[str, str]:
    """
    Globally align two sequences and return the aligned strings with gaps.

    Returns:
        Tuple of aligned sequence strings: (aligned_s1, aligned_s2)
    """
    if not s1 or not s2:
        return s1, s2
    aln = aligner.align(s1, s2)[0]
    coords = aln.coordinates
    a_coords, b_coords = coords[0], coords[1]

    out_a, out_b = [], []
    a_pos, b_pos = a_coords[0], b_coords[0]

    for i in range(1, len(a_coords)):
        a_next, b_next = a_coords[i], b_coords[i]
        a_span, b_span = a_next - a_pos, b_next - b_pos

        if a_span > 0 and b_span > 0:
            out_a.append(s1[a_pos:a_next])
            out_b.append(s2[b_pos:b_next])
        elif a_span > 0 and b_span == 0:
            out_a.append(s1[a_pos:a_next])
            out_b.append("-" * a_span)
        elif a_span == 0 and b_span > 0:
            out_a.append("-" * b_span)
            out_b.append(s2[b_pos:b_next])

        a_pos, b_pos = a_next, b_next

    return "".join(out_a), "".join(out_b)


def identity_ignore_gaps(aln_a: str, aln_b: str) -> Tuple[int, int, int, float]:
    """
    Calculate pairwise identity while ignoring gap-containing columns.

    Only positions where both aligned sequences have non-gap characters are
    compared.

    Returns:
        compared, matches, mismatches, identity_fraction
    """
    #### Compare only positions where BOTH are non-gap.
    #### Returns: compared, matches, mismatches, identity_fraction
    compared = matches = mismatches = 0
    for aa, bb in zip(aln_a, aln_b):
        if aa == "-" or bb == "-":
            continue
        compared += 1
        if aa == bb:
            matches += 1
        else:
            mismatches += 1

    identity = (matches / compared) if compared > 0 else 0.0
    return compared, matches, mismatches, identity


def out_csv_name(translated_filename: str) -> str:
    """
    Build the output CSV filename from a translated FASTA filename.
    """
    base = translated_filename.replace(".faa_protein.faa", "")
    return f"{base}.frameFAIL.csv"


def process_one_file(
    tpath: str,
    original_dir: str,
    out_dir: str,
    threshold: float,
    aligner: PairwiseAligner,
) -> str:
    """
    Process one translated FASTA file against its matching original FASTA file.

    Args:
        tpath: Path to the translated FASTA file.
        original_dir: Directory containing original FASTA files.
        out_dir: Directory where failure CSVs should be written.
        threshold: Minimum identity fraction required to pass.
        aligner: Configured PairwiseAligner object.

    Returns:
        Output CSV path if a failure CSV was written; otherwise "".

    Notes:
        A CSV is written only when at least one failure exists.
        Output rows include failures only; passing rows are not written.
    """
    #### Returns output path if written, else "".
    ####
    #### Writes a CSV ONLY when at least one "failure" exists. Failures include:
    ####   1) original fasta missing for this translated file
    ####   2) species_key present in translated but missing in original
    ####   3) pairwise identity < threshold (or compared==0)
    #### Output CSV contains ONLY failing rows (no YES rows).

    tbase = os.path.basename(tpath)
    obase = tbase.replace(".faa_protein.faa", "_ORIGINALSEQS.faa")
    opath = os.path.join(original_dir, obase)

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, out_csv_name(tbase))

    fieldnames = [
        "species_key",
        "translated_id",
        "original_id",
        "translated_len",
        "original_len",
        "nongap_compared",
        "matches",
        "mismatches",
        "identity_fraction",
        "identity_percent",
        "passes_threshold",
        "note",
    ]

    fail_rows: List[dict] = []

    #### Failure type 1: missing original fasta file
    if not os.path.exists(opath):
        fail_rows.append(
            {
                "species_key": "",
                "translated_id": "",
                "original_id": "",
                "translated_len": "",
                "original_len": "",
                "nongap_compared": 0,
                "matches": 0,
                "mismatches": 0,
                "identity_fraction": f"{0.0:.6f}",
                "identity_percent": f"{0.0:.2f}",
                "passes_threshold": "NO",
                "note": f"original fasta missing: {opath}",
            }
        )

        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            for row in fail_rows:
                w.writerow(row)

        print(f"{tbase}: FAIL (missing original) wrote={out_path}")
        return out_path

    trans = load_by_species_key(tpath)
    orig = load_by_species_key(opath)

    n_compared_pairs = 0
    n_missing_in_original = 0

    for skey, (tid, tseq) in trans.items():
        #### Failure type 2: species_key missing in original
        if skey not in orig:
            n_missing_in_original += 1
            fail_rows.append(
                {
                    "species_key": skey,
                    "translated_id": tid,
                    "original_id": "",
                    "translated_len": len(tseq),
                    "original_len": "",
                    "nongap_compared": 0,
                    "matches": 0,
                    "mismatches": 0,
                    "identity_fraction": f"{0.0:.6f}",
                    "identity_percent": f"{0.0:.2f}",
                    "passes_threshold": "NO",
                    "note": "species_key not present in original file",
                }
            )
            continue

        oid, oseq = orig[skey]
        aln_t, aln_o = align_to_gapped_strings(aligner, tseq, oseq)
        compared, matches, mismatches, identity = identity_ignore_gaps(aln_t, aln_o)
        n_compared_pairs += 1

        #### Failure type 3: identity below threshold, or compared==0
        passes = (identity >= threshold) if compared > 0 else False
        if not passes:
            fail_rows.append(
                {
                    "species_key": skey,
                    "translated_id": tid,
                    "original_id": oid,
                    "translated_len": len(tseq),
                    "original_len": len(oseq),
                    "nongap_compared": compared,
                    "matches": matches,
                    "mismatches": mismatches,
                    "identity_fraction": f"{identity:.6f}",
                    "identity_percent": f"{identity * 100.0:.2f}",
                    "passes_threshold": "NO",
                    "note": "",
                }
            )

    #### ONLY write if there is at least one failing row
    if len(fail_rows) == 0:
        return ""

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in fail_rows:
            w.writerow(row)

    print(
        f"{tbase}: compared_pairs={n_compared_pairs}, missing_in_original={n_missing_in_original}, "
        f"FAIL_ROWS={len(fail_rows)} wrote={out_path}"
    )
    return out_path


def main():
    """
    Parse command-line arguments and run the frame-failure comparison workflow.

    The script can process either:
        1) a single translated FASTA file using --infile, or
        2) all translated FASTA files in --translated_dir.
    """
    ap = argparse.ArgumentParser(
        #### no triple quotes; keep help simple
        description="Write a CSV only when failures exist; output contains only failing rows."
    )
    ap.add_argument("--translated_dir", required=True)
    ap.add_argument("--original_dir", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument(
        "--threshold",
        type=float,
        default=0.90,
        help="Identity threshold as a fraction (0.90 = 90% identity).",
    )

    ap.add_argument(
        "--infile",
        default=None,
        help="Single translated filename or full path (useful for SLURM array jobs).",
    )
    args = ap.parse_args()

    aligner = build_aligner()

    if args.infile:
        tpath = args.infile
        if not os.path.isabs(tpath):
            tpath = os.path.join(args.translated_dir, tpath)
        _ = process_one_file(
            tpath=tpath,
            original_dir=args.original_dir,
            out_dir=args.out_dir,
            threshold=args.threshold,
            aligner=aligner,
        )
        return

    translated_files = sorted(
        os.path.join(args.translated_dir, f)
        for f in os.listdir(args.translated_dir)
        if f.endswith(".faa_protein.faa")
    )

    for tpath in translated_files:
        _ = process_one_file(
            tpath=tpath,
            original_dir=args.original_dir,
            out_dir=args.out_dir,
            threshold=args.threshold,
            aligner=aligner,
        )


if __name__ == "__main__":
    main()
```
