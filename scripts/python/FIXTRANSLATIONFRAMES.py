```python
#!/usr/bin/env python3
"""
Frame fixer for translated protein FASTA files.

This script copies all translated protein FASTAs into a FRAME_FIXED-style output
directory. If a matching *.frameFAIL.csv exists, only the sequences flagged as
failures are replaced.

For each failing sequence, the script:
    1) Finds the matching CDS sequence.
    2) Finds the matching original protein sequence.
    3) Translates the CDS in all six reading frames.
    4) Compares each translated frame to the original protein.
    5) Replaces the failed translated sequence with the highest-identity frame.

Files with no frameFAIL CSV are copied unchanged.

Important:
    Failed sequences are replaced with the best available frame translation if
    one exists, even if the best identity does not improve or remains below the
    threshold. Warnings are printed in those cases.
"""

import argparse
import csv
import os
import re
import sys
import glob
import shutil
from typing import Dict, Tuple, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


def clean_seq(s: str) -> str:
    """
    Normalize a protein sequence.

    Removes whitespace, converts to uppercase, and removes trailing stop
    characters.
    """
    s = re.sub(r"\s+", "", str(s)).upper()
    #### drop terminal stop(s)
    while s.endswith("*"):
        s = s[:-1]
    return s


def species_key_before_first_underscore_from_id(seq_id: str) -> str:
    """
    Extract a species key from a FASTA record ID.

    Uses the first whitespace-delimited token, then keeps everything before the
    first underscore.
    """
    token = seq_id.split()[0]
    return token.split("_", 1)[0]


def build_aligner() -> PairwiseAligner:
    """
    Build the global pairwise aligner used for protein identity comparisons.
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
    Globally align two sequences and return aligned strings with gap characters.
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
    Calculate sequence identity while ignoring columns where either sequence has
    a gap.

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


def wrap60(s: str, width: int = 60) -> str:
    """
    Wrap a sequence string to fixed-width FASTA-style lines.
    """
    return "\n".join(s[i:i + width] for i in range(0, len(s), width))


def find_matching_framefail_csv(framefail_dir: str, translated_basename: str) -> Optional[str]:
    """
    Find the frameFAIL CSV associated with a translated protein FASTA.

    Matching is based on the translated filename prefix before
    .faa_protein.faa.
    """
    #### translated_basename is the filename only (not path)
    #### We match by the "start" prefix = basename without .faa_protein.faa
    prefix = translated_basename.replace(".faa_protein.faa", "")
    pattern = os.path.join(framefail_dir, prefix + "*.frameFAIL.csv")
    hits = sorted(glob.glob(pattern))
    if not hits:
        return None
    #### If multiple exist, take the first sorted; print warning
    if len(hits) > 1:
        print(f"WARNING: multiple frameFAIL CSVs match {pattern}; using {hits[0]}", file=sys.stderr)
    return hits[0]


def translated_to_cds_filename(translated_basename: str) -> str:
    """
    Convert a translated protein FASTA filename to the expected CDS filename.
    """
    #### 1000_... .faa_protein.faa  ->  1000_... .faa_CDS.fasta
    base = translated_basename.replace(".faa_protein.faa", "")
    return base + ".faa_CDS.fasta"


def translated_to_original_protein_filename(translated_basename: str) -> str:
    """
    Convert a translated protein FASTA filename to the expected original protein
    FASTA filename.
    """
    #### 1000_... .faa_protein.faa  ->  1000_... _ORIGINALSEQS.faa
    base = translated_basename.replace(".faa_protein.faa", "")
    return base + "_ORIGINALSEQS.faa"


def load_fasta_by_id(path: str) -> Dict[str, Seq]:
    """
    Load a FASTA file into a dictionary keyed by record ID.
    """
    d: Dict[str, Seq] = {}
    for rec in SeqIO.parse(path, "fasta"):
        d[rec.id] = rec.seq
    return d


def translate_six_frames(cds_seq: Seq, table: int = 1) -> List[Tuple[str, str]]:
    """
    Translate a CDS sequence in all six reading frames.

    Returns:
        List of frame_label, cleaned_protein_sequence pairs.

    Frames:
        +1, +2, +3, -1, -2, -3
    """
    #### returns list of (frame_label, protein_string_cleaned)
    #### frames: +1,+2,+3,-1,-2,-3
    frames: List[Tuple[str, str]] = []

    #### forward
    s = str(cds_seq).upper().replace(" ", "").replace("\n", "")
    s = re.sub(r"[^ACGTN]", "N", s)

    for off in (0, 1, 2):
        if len(s) - off < 3:
            prot = ""
        else:
            trimmed_len = (len(s) - off) - ((len(s) - off) % 3)
            sub = Seq(s[off:off + trimmed_len])
            prot = clean_seq(sub.translate(table=table, to_stop=False))
        frames.append((f"+{off+1}", prot))

    #### reverse complement
    rc = str(Seq(s).reverse_complement())
    for off in (0, 1, 2):
        if len(rc) - off < 3:
            prot = ""
        else:
            trimmed_len = (len(rc) - off) - ((len(rc) - off) % 3)
            sub = Seq(rc[off:off + trimmed_len])
            prot = clean_seq(sub.translate(table=table, to_stop=False))
        frames.append((f"-{off+1}", prot))

    return frames


def best_frame_translation(
    cds_seq: Seq,
    original_protein: str,
    aligner: PairwiseAligner,
    table: int,
) -> Tuple[str, str, float, str, str]:
    """
    Find the six-frame CDS translation that best matches the original protein.

    Each translated frame is globally aligned to the original protein, then
    scored using gap-ignored identity.

    Returns:
        best_frame_label, best_protein, best_identity, best_aln_a, best_aln_b
    """
    #### returns (best_frame_label, best_protein, best_identity, best_aln_a, best_aln_b)
    best_frame = "NA"
    best_prot = ""
    best_ident = -1.0
    best_aln_a = ""
    best_aln_b = ""

    for frame_label, prot in translate_six_frames(cds_seq, table=table):
        if not prot:
            continue
        aln_a, aln_b = align_to_gapped_strings(aligner, prot, original_protein)
        _, _, _, ident = identity_ignore_gaps(aln_a, aln_b)
        if ident > best_ident:
            best_ident = ident
            best_frame = frame_label
            best_prot = prot
            best_aln_a = aln_a
            best_aln_b = aln_b

    if best_ident < 0:
        #### no usable translation
        best_ident = 0.0

    return best_frame, best_prot, best_ident, best_aln_a, best_aln_b


def current_identity(
    current_translated: str,
    original_protein: str,
    aligner: PairwiseAligner,
) -> Tuple[float, str, str]:
    """
    Compare the current translated protein sequence to the original protein.

    Returns:
        identity_fraction, aligned_current, aligned_original
    """
    aln_a, aln_b = align_to_gapped_strings(aligner, current_translated, original_protein)
    _, _, _, ident = identity_ignore_gaps(aln_a, aln_b)
    return ident, aln_a, aln_b


def parse_framefail_csv(csv_path: str) -> List[dict]:
    """
    Read a frameFAIL CSV and return rows where passes_threshold is NO.
    """
    #### returns list of failing rows (already should be only fails, but we enforce NO)
    fails: List[dict] = []
    with open(csv_path, "r", newline="") as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            if row.get("passes_threshold", "").strip().upper() == "NO":
                fails.append(row)
    return fails


def main():
    """
    Run the frame-fixing workflow.

    Processes either one translated protein FASTA supplied with --infile, or all
    *.faa_protein.faa files in --translated_dir.

    For files without a matching frameFAIL CSV, the input FASTA is copied
    unchanged. For files with failures, only failed sequences are replaced using
    the best six-frame CDS translation.
    """
    ap = argparse.ArgumentParser(
        description="Frame fixer: copy all translated protein FASTAs into FRAME_FIXED, and replace only sequences flagged in *.frameFAIL.csv by best of 6-frame CDS translations."
    )
    ap.add_argument("--translated_dir", required=True)
    ap.add_argument("--framefail_dir", required=True)
    ap.add_argument("--cds_dir", required=True)
    ap.add_argument("--original_protein_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--threshold", type=float, default=0.90)
    ap.add_argument("--genetic_code_table", type=int, default=1)
    ap.add_argument("--infile", default=None, help="Single translated protein fasta (basename or full path). If provided, only process this one file.")

    args = ap.parse_args()

    translated_dir = args.translated_dir
    framefail_dir = args.framefail_dir
    cds_dir = args.cds_dir
    original_protein_dir = args.original_protein_dir
    out_dir = args.out_dir
    threshold = args.threshold
    table = args.genetic_code_table

    os.makedirs(out_dir, exist_ok=True)

    aligner = build_aligner()

    #### choose which files to process
    if args.infile:
        #### allow basename or full path
        if os.path.isabs(args.infile):
            fname = os.path.basename(args.infile)
            translated_files = [fname]
        else:
            translated_files = [args.infile]
    else:
        translated_files = sorted(
            f for f in os.listdir(translated_dir)
            if f.endswith(".faa_protein.faa")
        )

    if not translated_files:
        print(f"ERROR: no *.faa_protein.faa files found in {translated_dir}", file=sys.stderr)
        sys.exit(2)

    n_total = 0
    n_copied_as_is = 0
    n_fixed_files = 0
    n_total_replaced = 0
    n_warn = 0

    for fname in translated_files:
        n_total += 1
        in_faa = os.path.join(translated_dir, fname)

        if not os.path.exists(in_faa):
            print(f"ERROR: translated protein file not found: {in_faa}", file=sys.stderr)
            n_warn += 1
            continue

        out_faa = os.path.join(out_dir, fname)

        csv_path = find_matching_framefail_csv(framefail_dir, fname)

        #### If no frameFAIL csv: copy as-is
        if csv_path is None:
            shutil.copy2(in_faa, out_faa)
            n_copied_as_is += 1
            continue

        #### We have failures: load failing rows
        fail_rows = parse_framefail_csv(csv_path)

        #### If somehow csv exists but has no NO rows, treat as clean copy
        if len(fail_rows) == 0:
            shutil.copy2(in_faa, out_faa)
            n_copied_as_is += 1
            continue

        #### Locate matching CDS fasta and original protein fasta
        cds_file = os.path.join(cds_dir, translated_to_cds_filename(fname))
        orig_prot_file = os.path.join(original_protein_dir, translated_to_original_protein_filename(fname))

        if not os.path.exists(cds_file):
            print(f"WARNING: missing CDS file for {fname}: {cds_file}  -> copying as-is", file=sys.stderr)
            shutil.copy2(in_faa, out_faa)
            n_warn += 1
            continue

        if not os.path.exists(orig_prot_file):
            print(f"WARNING: missing original protein file for {fname}: {orig_prot_file}  -> copying as-is", file=sys.stderr)
            shutil.copy2(in_faa, out_faa)
            n_warn += 1
            continue

        #### Load translated proteins (as SeqRecords to preserve order and IDs)
        translated_records = list(SeqIO.parse(in_faa, "fasta"))
        tr_by_id: Dict[str, int] = {rec.id: i for i, rec in enumerate(translated_records)}

        #### Load CDS records and original protein records keyed by exact rec.id
        cds_by_id = load_fasta_by_id(cds_file)
        orig_by_id = load_fasta_by_id(orig_prot_file)

        replaced_in_this_file = 0

        #### For each failing row, fix that translated_id if possible
        for row in fail_rows:
            t_id = row.get("translated_id", "").strip()
            species_key = row.get("species_key", "").strip()

            if not t_id:
                print(f"WARNING: {fname}: CSV row missing translated_id; skipping row", file=sys.stderr)
                n_warn += 1
                continue

            #### The original protein header should match by the same ID in your examples.
            #### If not found, fallback to species_key prefix matching.
            orig_id = row.get("original_id", "").strip()

            #### Find translated record
            if t_id not in tr_by_id:
                #### fallback: try species_key prefix match among translated
                found = None
                for rec_id in tr_by_id.keys():
                    if species_key and species_key_before_first_underscore_from_id(rec_id) == species_key:
                        found = rec_id
                        break
                if found is None:
                    print(f"WARNING: {fname}: cannot find translated record for {t_id} (species_key={species_key}); skipping", file=sys.stderr)
                    n_warn += 1
                    continue
                t_id = found

            tr_idx = tr_by_id[t_id]
            current_tr = clean_seq(translated_records[tr_idx].seq)

            #### Find CDS record
            cds_seq = None
            if t_id in cds_by_id:
                cds_seq = cds_by_id[t_id]
            else:
                #### fallback: species_key prefix match in CDS
                for cid, cseq in cds_by_id.items():
                    if species_key and species_key_before_first_underscore_from_id(cid) == species_key:
                        cds_seq = cseq
                        break

            if cds_seq is None:
                print(f"WARNING: {fname}: cannot find CDS for {t_id} (species_key={species_key}); leaving sequence unchanged", file=sys.stderr)
                n_warn += 1
                continue

            #### Find original protein sequence
            original_prot_seq = None

            #### preferred: use original_id from CSV if present
            if orig_id and orig_id in orig_by_id:
                original_prot_seq = clean_seq(orig_by_id[orig_id])
            elif t_id in orig_by_id:
                original_prot_seq = clean_seq(orig_by_id[t_id])
            else:
                #### fallback: species_key prefix match in originals
                for oid2, oseq2 in orig_by_id.items():
                    if species_key and species_key_before_first_underscore_from_id(oid2) == species_key:
                        original_prot_seq = clean_seq(oseq2)
                        break

            if original_prot_seq is None:
                print(f"WARNING: {fname}: cannot find original protein for {t_id} (species_key={species_key}); leaving unchanged", file=sys.stderr)
                n_warn += 1
                continue

            #### Compute current identity
            cur_ident, cur_aln_a, cur_aln_b = current_identity(current_tr, original_prot_seq, aligner)

            #### Find best frame translation
            best_frame, best_prot, best_ident, best_aln_a, best_aln_b = best_frame_translation(
                cds_seq=cds_seq,
                original_protein=original_prot_seq,
                aligner=aligner,
                table=table,
            )

            #### Decide replacement: you requested highest identity, but warn if no improvement
            if best_ident <= cur_ident + 1e-12:
                print(
                    f"WARNING: {fname}: {t_id}: no frame improved identity (current={cur_ident*100:.2f}%, best={best_ident*100:.2f}% frame={best_frame})",
                    file=sys.stderr
                )
                n_warn += 1

            #### Replace regardless with best (highest identity) if we have one
            if best_prot:
                translated_records[tr_idx].seq = Seq(best_prot)
                replaced_in_this_file += 1
                n_total_replaced += 1

            #### Warn and print best alignment if still below threshold
            if best_ident < threshold - 1e-12:
                n_warn += 1
                print(
                    f"WARNING: {fname}: {t_id}: best frame still below threshold {threshold*100:.1f}% "
                    f"(best={best_ident*100:.2f}% frame={best_frame})",
                    file=sys.stderr
                )
                print("#### BEST ALIGNMENT (candidate vs original) ####", file=sys.stderr)
                print(f"#### file={fname}", file=sys.stderr)
                print(f"#### id={t_id}", file=sys.stderr)
                print(f"#### frame={best_frame}", file=sys.stderr)
                print(f"#### identity={best_ident*100:.2f}%", file=sys.stderr)
                print(wrap60(best_aln_a), file=sys.stderr)
                print(wrap60(best_aln_b), file=sys.stderr)
                print("#### END ALIGNMENT ####", file=sys.stderr)

        #### Write output file (fixed or partially fixed)
        #### You requested FRAME_FIXED contain ALL files, so we always write it here
        SeqIO.write(translated_records, out_faa, "fasta")
        n_fixed_files += 1
        print(f"{fname}: wrote {out_faa}  replaced={replaced_in_this_file}  fails_in_csv={len(fail_rows)}")

    print(
        f"DONE. total_files={n_total} copied_as_is={n_copied_as_is} processed_with_csv={n_fixed_files} "
        f"total_sequences_replaced={n_total_replaced} warnings={n_warn}"
    )

if __name__ == "__main__":
    main()
```
