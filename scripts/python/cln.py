```python
#!/usr/bin/env python3
"""
cln.py

Codon-only cleaner for HyPhy prep.

Input:
    A directory containing filelist.txt, with one FASTA path per line.
    Relative paths are allowed and are resolved relative to --indir.

For each FASTA, this script:
    1) Writes a CSV backup of full, unedited headers into --backupdir.
    2) Writes a cleaned FASTA into --outdir.
    3) Truncates FASTA headers at the first delimiter, default "_".
    4) Replaces stop codons TAA/TAG/TGA with gaps "---".

Designed for SLURM arrays:
    sbatch --array=1-N ...
    Uses SLURM_ARRAY_TASK_ID, 1-based, to select the Nth entry in filelist.txt.
"""

## cln.py
## Codon-only cleaner for HyPhy prep.
## Input: a directory that contains filelist.txt (one FASTA path per line; relative paths allowed)
## Does, per FASTA:
##   1) writes a CSV (one row per FULL unedited header) into --backupdir
##   2) writes a NEW FASTA into --outdir with headers truncated at first "_" and filename suffix "_cln"
##   3) replaces STOP CODONS (TAA/TAG/TGA) with gaps '---' (CODON MODE ONLY; never protein)
##
## Designed for SLURM arrays:
##   sbatch --array=1-N ...
##   uses $SLURM_ARRAY_TASK_ID (1-based) to pick the Nth entry in filelist.txt

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Iterator, Tuple, List


def iter_fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Iterate through FASTA records from a file.

    Args:
        path: Path to the input FASTA file.

    Yields:
        Tuples of (full_header_without_>, sequence_concatenated).

    Notes:
        Blank lines are ignored. Sequence lines are concatenated into one string
        per record.
    """
    ## Yields (full_header_without_>, sequence_concatenated)
    header = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]  ## full header (unedited)
                seq_chunks = []
            else:
                seq_chunks.append(line.replace(" ", "").replace("\t", ""))
        if header is not None:
            yield header, "".join(seq_chunks)


def wrap(seq: str, width: int) -> str:
    """
    Wrap a sequence string to a fixed FASTA line width.

    Args:
        seq: Sequence string to wrap.
        width: Number of characters per output line.

    Returns:
        Wrapped sequence string.
    """
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def clean_id(full_header: str, delim: str) -> str:
    """
    Create a HyPhy-safe taxon ID from a full FASTA header.

    Uses the first whitespace-delimited token, then truncates that token at the
    first delimiter occurrence.

    Args:
        full_header: Full FASTA header without the leading ">".
        delim: Delimiter used to truncate IDs.

    Returns:
        Cleaned sequence/taxon ID.
    """
    ## HyPhy name = first token up to whitespace, then truncate at first delim
    first_token = full_header.split()[0]
    if delim in first_token:
        return first_token.split(delim, 1)[0]
    return first_token


def ensure_unique(base: str, used: set[str]) -> str:
    """
    Ensure a cleaned taxon label is unique within one output FASTA.

    If the base label has already been used, appends _2, _3, etc.

    Args:
        base: Proposed cleaned label.
        used: Set of labels already used in this output FASTA.

    Returns:
        Unique label.
    """
    ## Guarantee unique taxa labels within an output FASTA
    if base not in used:
        used.add(base)
        return base
    k = 2
    while f"{base}_{k}" in used:
        k += 1
    out = f"{base}_{k}"
    used.add(out)
    return out


def read_filelist(list_path: Path, base_dir: Path) -> List[Path]:
    """
    Read FASTA targets from a file list.

    Args:
        list_path: Path to filelist.txt.
        base_dir: Directory used to resolve relative paths.

    Returns:
        List of resolved FASTA paths.

    Notes:
        Blank lines and comment lines beginning with "#" are ignored.
    """
    ## Reads one file path per line, ignores blank lines and comments starting with '#'
    ## If a line is relative, it is resolved relative to base_dir (the input directory).
    items: List[Path] = []
    with list_path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            p = Path(line)
            if not p.is_absolute():
                p = (base_dir / p).resolve()
            items.append(p)
    return items


def is_codonish(seq: str) -> bool:
    """
    Check whether a sequence looks like codon/nucleotide data.

    Allows nucleotide bases, common ambiguity codes, gaps, missing data
    characters, and "!" characters. Rejects "*" because this script is not
    intended for protein FASTA files.

    Args:
        seq: Sequence string to check.

    Returns:
        True if the sequence has an allowed codon-mode alphabet; otherwise False.
    """
    ## Allow nucleotide letters + common ambiguity + gaps/missing. No '*' allowed in codon mode.
    ## If '*' exists, we fail loudly (you said never proteins).
    if "*" in seq:
        return False
    return re.fullmatch(r"[ACGTUNRYWSKMBDHVacgtunrywskmbdhv\-\?\.\!]*", seq) is not None


def replace_stop_codons_with_gaps_codon_only(seq: str) -> str:
    """
    Replace in-frame stop codons with gap codons.

    Args:
        seq: Codon-aligned nucleotide sequence.

    Returns:
        Uppercase sequence with TAA, TAG, and TGA replaced by "---".

    Raises:
        SystemExit: If the sequence contains non-codon characters, contains "*",
        or has a length that is not divisible by 3.

    Notes:
        Stop codons interrupted by gaps are not detected.
    """
    ## Codon-only:
    ## - requires nucleotide-ish alphabet
    ## - requires length % 3 == 0
    ## - replaces TAA/TAG/TGA with '---'
    ## Notes:
    ## - Works in uppercase for matching; outputs uppercase.
    ## - If your alignments include lots of gaps, a stop codon spanning gaps will NOT be detected.
    if not is_codonish(seq):
        raise SystemExit("ERROR: sequence contains non-codon characters (or '*'). This script is codon-only.")

    if len(seq) % 3 != 0:
        raise SystemExit(f"ERROR: sequence length not divisible by 3 in codon mode (len={len(seq)}).")

    s = seq.upper()
    stops = {"TAA", "TAG", "TGA"}
    out = []
    for i in range(0, len(s), 3):
        codon = s[i : i + 3]
        if codon in stops:
            out.append("---")
        else:
            out.append(codon)
    return "".join(out)


def process_one(fasta_path: Path, outdir: Path, backupdir: Path, suffix: str, delim: str, width: int) -> None:
    """
    Clean one FASTA file and write outputs.

    Args:
        fasta_path: Input FASTA path.
        outdir: Directory for cleaned FASTA output.
        backupdir: Directory for header-backup CSV.
        suffix: Suffix appended to the output FASTA filename stem.
        delim: Delimiter used to truncate FASTA IDs.
        width: FASTA line wrapping width.

    Outputs:
        - Header backup CSV in backupdir.
        - Cleaned FASTA in outdir.

    Raises:
        SystemExit: If the input FASTA does not exist or contains no records.
    """
    if not fasta_path.exists():
        raise SystemExit(f"ERROR: input FASTA not found: {fasta_path}")

    records = list(iter_fasta_records(fasta_path))
    if not records:
        raise SystemExit(f"ERROR: no FASTA records found in: {fasta_path}")

    outdir.mkdir(parents=True, exist_ok=True)
    backupdir.mkdir(parents=True, exist_ok=True)

    ## 1) CSV backup with FULL headers
    csv_path = backupdir / f"{fasta_path.stem}_headers.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as out_csv:
        w = csv.writer(out_csv)
        w.writerow(["record_index", "original_full_header"])
        for i, (hdr, _seq) in enumerate(records, start=1):
            w.writerow([i, hdr])

    ## 2 + 3) output FASTA with cleaned IDs and stop codons -> gaps
    out_fasta = outdir / f"{fasta_path.stem}{suffix}{fasta_path.suffix}"
    used: set[str] = set()

    with out_fasta.open("w", encoding="utf-8") as out_fh:
        for hdr, seq in records:
            cid = ensure_unique(clean_id(hdr, delim=delim), used)
            seq2 = replace_stop_codons_with_gaps_codon_only(seq)
            out_fh.write(f">{cid}\n")
            out_fh.write(wrap(seq2, width=width))
            out_fh.write("\n")

    print(f"OK\tIN={fasta_path}\tOUT={out_fasta}\tCSV={csv_path}")


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed argparse Namespace.
    """
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument(
        "--indir",
        type=Path,
        required=True,
        help="Input directory containing filelist.txt and (typically) the FASTA files",
    )
    ap.add_argument("--outdir", type=Path, required=True, help="Directory for cleaned FASTA outputs")
    ap.add_argument("--backupdir", type=Path, required=True, help="Directory for per-file header CSV backups")

    ap.add_argument("--filelist-name", default="filelist.txt", help='Name of file list inside --indir (default: "filelist.txt")')

    ap.add_argument("--suffix", default="_cln", help='Suffix appended to output filename stem (default: "_cln")')
    ap.add_argument("--delim", default="_", help='Delimiter to truncate IDs (default: "_")')
    ap.add_argument("--line-width", type=int, default=60, help="FASTA wrap width (default: 60)")

    ap.add_argument(
        "--array",
        action="store_true",
        help="If set, process exactly one file determined by SLURM_ARRAY_TASK_ID (1-based) from filelist.txt",
    )

    return ap.parse_args()


def main() -> None:
    """
    Run the codon-only FASTA cleaning workflow.

    Reads filelist.txt from --indir, optionally selects one target using
    SLURM_ARRAY_TASK_ID, and processes each selected FASTA.
    """
    args = parse_args()

    indir: Path = args.indir
    if not indir.is_dir():
        raise SystemExit(f"ERROR: --indir is not a directory: {indir}")

    list_path = indir / args.filelist_name
    if not list_path.exists():
        raise SystemExit(f"ERROR: expected file list not found: {list_path}")

    targets = read_filelist(list_path, base_dir=indir)
    if not targets:
        raise SystemExit(f"ERROR: no targets found in {list_path}")

    if args.array:
        tid = os.environ.get("SLURM_ARRAY_TASK_ID")
        if tid is None:
            raise SystemExit("ERROR: --array set but SLURM_ARRAY_TASK_ID is not defined in environment")
        try:
            idx = int(tid)
        except ValueError:
            raise SystemExit(f"ERROR: SLURM_ARRAY_TASK_ID is not an integer: {tid}")

        if idx < 1 or idx > len(targets):
            raise SystemExit(f"ERROR: SLURM_ARRAY_TASK_ID={idx} out of range (1..{len(targets)})")

        targets = [targets[idx - 1]]

    for fp in targets:
        process_one(
            fasta_path=fp,
            outdir=args.outdir,
            backupdir=args.backupdir,
            suffix=args.suffix,
            delim=args.delim,
            width=args.line_width,
        )


if __name__ == "__main__":
    main()
```
