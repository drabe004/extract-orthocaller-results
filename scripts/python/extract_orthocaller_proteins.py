```python
#!/usr/bin/env python3
"""
Extract Orthocaller-selected protein sequences from GeneRax alignment FASTAs.

This script starts from a master Orthocaller summary, identifies the selected
orthogroup entries, maps those entries back to the corresponding GeneRax
alignment FASTA files, and extracts the matching protein sequences.

Inputs:
    --master_summary
        File listing target orthogroups. Each useful line is expected to begin
        with an orthogroup ID before the first ":".

    --orthocaller_base
        Base directory containing Orthocaller output subdirectories.

    --generax_key
        CSV-like key mapping GeneRax orthogroup IDs to alignment FASTA names.

    --aln_dir
        Directory containing GeneRax alignment FASTA files.

Outputs:
    - extraction_report.tsv
    - Either one combined FASTA or one extracted FASTA per orthogroup.

Typical use:
    Recover the exact protein sequences selected by Orthocaller for downstream
    annotation, QC, or comparison against original GeneRax alignments.
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple


def read_master_summary(master_summary_path: Path) -> List[str]:
    """
    Read target orthogroup IDs from a master summary file.

    Each non-empty, non-comment line is parsed by taking the text before the
    first colon.

    Args:
        master_summary_path: Path to the master summary file.

    Returns:
        List of orthogroup IDs.
    """
    ogs = []

    with master_summary_path.open() as fh:
        for line in fh:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            og = line.split(":", 1)[0].strip()

            if og:
                ogs.append(og)

    return ogs


def read_generax_key(key_path: Path) -> Dict[str, str]:
    """
    Read the GeneRax alignment key file.

    Expected format:
        fasta_filename,generax_orthogroup_id

    Args:
        key_path: Path to the GeneRax key file.

    Returns:
        Dictionary mapping GeneRax orthogroup ID -> FASTA filename.
    """
    mapping = {}

    with key_path.open() as fh:
        for raw in fh:
            line = raw.strip()

            if not line or line.startswith("#"):
                continue

            parts = [p.strip() for p in line.split(",")]

            if len(parts) != 2:
                continue

            fname, generax = parts
            mapping[generax] = fname

    return mapping


def og_dir_from_og_full(og_full: str) -> str:
    """
    Convert a full Orthocaller gene ID into its parent orthogroup directory name.

    Args:
        og_full: Full orthogroup/gene identifier.

    Returns:
        Parent orthogroup directory name.
    """
    if "-Gene-" in og_full:
        return og_full.split("-Gene-", 1)[0]
    return og_full


def csv_row_candidates(og_full: str) -> List[str]:
    """
    Generate possible Orthocaller CSV row IDs for a target orthogroup.

    Master summary/output may use 107b_generax-Gene-18, but the internal
    orthogroup CSV row may still use 107_generax-Gene-18. This function tries
    both while preserving og_full for output.

    Args:
        og_full: Full orthogroup/gene ID from the master summary.

    Returns:
        Ordered list of candidate row IDs to try.
    """
    """
    Master summary/output may use 107b_generax-Gene-18,
    but the internal orthogroup CSV row may still use 107_generax-Gene-18.
    Try both, while preserving og_full for output.
    """
    candidates = [og_full]

    if "b_generax" in og_full:
        candidates.append(og_full.replace("b_generax", "_generax"))

    return list(dict.fromkeys(candidates))


def read_orthocaller_row(orthocaller_csv: Path, target_og: str) -> List[str]:
    """
    Read selected sequence entries for one orthogroup from an Orthocaller CSV.

    Args:
        orthocaller_csv: Path to the Orthocaller orthogroup CSV.
        target_og: Orthogroup row ID to find.

    Returns:
        List of selected entry identifiers from the matching row.

    Raises:
        FileNotFoundError: If the Orthocaller CSV does not exist.
        ValueError: If the target orthogroup row is not found.
    """
    if not orthocaller_csv.exists():
        raise FileNotFoundError(f"Missing orthocaller csv: {orthocaller_csv}")

    with orthocaller_csv.open() as fh:
        reader = csv.reader(fh)

        for row in reader:
            if not row:
                continue

            og_id = row[0].strip()

            if og_id == target_og:
                payload = ",".join(row[1:]).strip()
                entries = re.split(r"\s+", payload) if payload else []
                return [e for e in entries if e]

    raise ValueError(f"Orthogroup {target_og} not found in {orthocaller_csv}")


def read_orthocaller_row_with_fallback(orthocaller_csv: Path, og_full: str) -> Tuple[List[str], str]:
    """
    Read an Orthocaller row, trying fallback row-name candidates if needed.

    Args:
        orthocaller_csv: Path to Orthocaller orthogroup CSV.
        og_full: Full orthogroup/gene ID from master summary.

    Returns:
        Tuple of:
            - selected entry identifiers
            - row ID that successfully matched

    Raises:
        Exception: Re-raises the last error if all candidates fail.
    """
    last_error = None

    for candidate in csv_row_candidates(og_full):
        try:
            entries = read_orthocaller_row(orthocaller_csv, candidate)
            return entries, candidate
        except Exception as e:
            last_error = e

    raise last_error


def parse_key_from_entry(entry: str) -> str:
    """
    Extract the searchable sequence key from an Orthocaller entry.

    If the entry contains underscores, the portion after the final underscore is
    used. Otherwise the full entry is returned.

    Args:
        entry: Orthocaller selected sequence entry.

    Returns:
        Sequence key used to search FASTA headers.
    """
    if "_" not in entry:
        return entry

    return entry.rsplit("_", 1)[-1]


def fasta_iter(fasta_path: Path):
    """
    Iterate through FASTA records.

    Args:
        fasta_path: Path to FASTA file.

    Yields:
        Tuples of (header_without_>, sequence).
    """
    header = None
    chunks = []

    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)

                header = line.strip()[1:]
                chunks = []
            else:
                chunks.append(line.strip())

    if header is not None:
        yield header, "".join(chunks)


def collect_matches(fasta_path: Path, target_keys: List[str]) -> Tuple[List[Tuple[str, str]], Dict[str, bool]]:
    """
    Collect FASTA records whose headers contain target sequence keys.

    Args:
        fasta_path: Path to GeneRax alignment FASTA.
        target_keys: Sequence keys parsed from Orthocaller entries.

    Returns:
        Tuple of:
            - list of matching (header, sequence) records
            - dictionary tracking whether each target key was seen
    """
    seen = {k: False for k in target_keys}
    matches = []

    for hdr, seq in fasta_iter(fasta_path):
        for k in target_keys:
            if k in hdr:
                matches.append((hdr, seq))
                seen[k] = True

    return matches, seen


def write_fasta_record(handle, header: str, seq: str):
    """
    Write one FASTA record to an open file handle.

    Args:
        handle: Open writable file handle.
        header: FASTA header without leading ">".
        seq: Sequence string.
    """
    handle.write(f">{header}\n")

    for i in range(0, len(seq), 60):
        handle.write(seq[i:i + 60] + "\n")


def main():
    """
    Run the Orthocaller-to-GeneRax sequence extraction workflow.

    For each target orthogroup:
        1) Locate the Orthocaller orthogroup CSV.
        2) Read selected sequence entries.
        3) Find the corresponding GeneRax alignment FASTA.
        4) Extract matching FASTA records.
        5) Write per-OG FASTAs or one combined FASTA.
        6) Write an extraction report.

    If --strict is set, missing files or missing sequence keys stop the run.
    Otherwise, errors are recorded in the report and processing continues.
    """
    ap = argparse.ArgumentParser(
        description="Extract Orthocaller-selected protein sequences from GeneRax alignment FASTAs."
    )

    ap.add_argument("--master_summary", required=True, type=Path)
    ap.add_argument("--orthocaller_base", required=True, type=Path)
    ap.add_argument("--generax_key", required=True, type=Path)
    ap.add_argument("--aln_dir", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--combined", action="store_true")
    ap.add_argument("--strict", action="store_true")

    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    report_path = args.out_dir / "extraction_report.tsv"
    combined_path = args.out_dir / "combined_extracted.faa"

    target_ogs = read_master_summary(args.master_summary)

    if not target_ogs:
        raise SystemExit("No OGs parsed from master summary.")

    grx_map = read_generax_key(args.generax_key)

    report_rows = []
    total_found = 0
    total_expected = 0

    combined_fh = combined_path.open("w") if args.combined else None

    try:
        for og_full in target_ogs:
            og_dir = og_dir_from_og_full(og_full)

            orthocaller_csv = args.orthocaller_base / og_dir / f"{og_dir}_orthogroup.csv"

            try:
                entries, matched_csv_row = read_orthocaller_row_with_fallback(orthocaller_csv, og_full)
            except Exception as e:
                report_rows.append([og_full, "0", "0", "NA", f"Error reading orthocaller row: {e}"])

                if args.strict:
                    raise

                continue

            if not entries:
                report_rows.append([og_full, "0", "0", "NA", "No entries in Orthocaller row"])

                if args.strict:
                    raise SystemExit(f"No entries found for {og_full}")

                continue

            keys = [parse_key_from_entry(e) for e in entries]
            total_expected += len(keys)

            if og_dir not in grx_map:
                report_rows.append([og_full, str(len(keys)), "0", "NA", f"Missing {og_dir} in GeneRax key"])

                if args.strict:
                    raise SystemExit(f"{og_dir} not found in GeneRax key")

                continue

            fasta_name = grx_map[og_dir]
            fasta_path = args.aln_dir / fasta_name

            if not fasta_path.exists():
                report_rows.append([og_full, str(len(keys)), "0", str(fasta_path), "FASTA missing"])

                if args.strict:
                    raise SystemExit(f"FASTA not found: {fasta_path}")

                continue

            matches, seen = collect_matches(fasta_path, keys)
            found = len({h for h, _ in matches})
            total_found += found

            if args.combined:
                for hdr, seq in matches:
                    write_fasta_record(combined_fh, f"{og_full}|{hdr}", seq)
            else:
                aln_stem = fasta_path.stem
                out_faa = args.out_dir / f"{og_full}__{aln_stem}.faa"

                with out_faa.open("w") as oh:
                    for hdr, seq in matches:
                        write_fasta_record(oh, hdr, seq)

            missing_keys = [k for k, ok in seen.items() if not ok]
            note = ""

            if matched_csv_row != og_full:
                note = f"CSV row matched by fallback: {matched_csv_row}"

            if missing_keys:
                missing_note = f"Missing {len(missing_keys)} keys; example: {missing_keys[0]}"
                note = f"{note}; {missing_note}" if note else missing_note

                if args.strict:
                    raise SystemExit(f"{og_full}: {len(missing_keys)} keys not found in {fasta_path.name}")

            report_rows.append([og_full, str(len(keys)), str(found), str(fasta_path), note])

    finally:
        if combined_fh:
            combined_fh.close()

    with report_path.open("w") as rh:
        rh.write("orthogroup\texpected_count\tfound_count\tfasta_path\tnote\n")

        for r in report_rows:
            rh.write("\t".join(r) + "\n")

    print(f"[OK] Wrote report: {report_path}")

    if args.combined:
        print(f"[OK] Wrote combined FASTA: {combined_path}")
    else:
        print(f"[OK] Wrote per-OG FASTAs under: {args.out_dir}")

    print(f"[SUMMARY] Found {total_found} / {total_expected} sequences across {len(target_ogs)} OGs.")


if __name__ == "__main__":
    main()
```