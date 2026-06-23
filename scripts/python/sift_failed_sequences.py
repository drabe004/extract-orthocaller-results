#!/usr/bin/env python3
"""
Remove failed sequences listed in frameFAIL CSV files from matching CDS and
protein FASTA files.

The script reads one or more *.frameFAIL.csv files, extracts failed translated
and original sequence IDs, removes matching records from protein and CDS FASTAs,
and writes sifted FASTA outputs to new directories.
"""

import argparse
import csv
import shutil
from pathlib import Path
from typing import Set, Tuple, Optional, List

from Bio import SeqIO

FASTA_EXTS_DEFAULT = [".fa", ".faa", ".fasta", ".fas", ".fna"]


def parse_args():
    """
    Parse command-line arguments for frame-failure sequence removal.
    """
    p = argparse.ArgumentParser(
        description="Remove failed sequences (from frameFAIL CSVs) from matching CDS and protein FASTAs, writing to new dirs."
    )
    p.add_argument("--fail_csv_dir", required=True, type=Path, help="Directory containing *.frameFAIL.csv files")
    p.add_argument("--cds_dir", required=True, type=Path, help="Directory containing CDS FASTAs")
    p.add_argument("--protein_dir", required=True, type=Path, help="Directory containing protein (translated CDS) FASTAs")
    p.add_argument("--out_cds_dir", required=True, type=Path, help="Output directory for sifted CDS FASTAs")
    p.add_argument("--out_protein_dir", required=True, type=Path, help="Output directory for sifted protein FASTAs")
    p.add_argument("--fasta_exts", nargs="+", default=FASTA_EXTS_DEFAULT,
                   help="FASTA extensions to consider (default: .fa .faa .fasta .fas .fna)")
    p.add_argument("--pattern", default="*.frameFAIL.csv", help="Glob for fail CSVs (default: *.frameFAIL.csv)")
    p.add_argument("--remove_if_passes_is", default="NO",
                   help="Only remove rows where passes_threshold == this value (default: NO). "
                        "If CSVs are fail-only, keep default.")
    p.add_argument("--dry_run", action="store_true", help="Do not write outputs; print what would happen.")

    # ---- ARRAY/SELECTOR OPTIONS ----
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--csv_index", type=int,
                      help="1-based index into the CSV list (use with SLURM_ARRAY_TASK_ID). Processes ONE CSV.")
    mode.add_argument("--single_csv", type=Path,
                      help="Path to a single .frameFAIL.csv to process (processes ONE CSV).")
    p.add_argument("--sort_csvs", action="store_true",
                   help="Sort CSV list deterministically before indexing (recommended for arrays).")

    return p.parse_args()


def record_matches_id(record, target_id: str) -> bool:
    """
    Check whether a FASTA record matches a target sequence ID.
    """
    # match record.id, record.name, or substring of description
    if record.id == target_id:
        return True
    if getattr(record, "name", "") == target_id:
        return True
    desc = getattr(record, "description", "")
    return target_id in desc


def find_matching_fasta(base_stem: str, search_dir: Path, fasta_exts: List[str]) -> Optional[Path]:
    """
    Find the FASTA file that corresponds to a frameFAIL CSV basename.
    """
    # exact stem + ext
    for ext in fasta_exts:
        cand = search_dir / f"{base_stem}{ext}"
        if cand.exists():
            return cand
    # prefix matches
    for ext in fasta_exts:
        hits = sorted(search_dir.glob(f"{base_stem}*{ext}"))
        if hits:
            return hits[0]
    return None


def load_fail_ids(fail_csv: Path, remove_if_passes_is: str) -> Tuple[Set[str], Set[str]]:
    """
    Load translated and original sequence IDs that should be removed.
    """
    translated_remove: Set[str] = set()
    original_remove: Set[str] = set()
    target = remove_if_passes_is.strip().upper()

    with fail_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            passes = (row.get("passes_threshold") or "").strip().upper()
            # If passes_threshold exists, filter; if missing/blank, treat as fail row
            if passes and passes != target:
                continue

            tid = (row.get("translated_id") or "").strip()
            oid = (row.get("original_id") or "").strip()
            if tid:
                translated_remove.add(tid)
            if oid:
                original_remove.add(oid)

    return translated_remove, original_remove


def sift_fasta(in_fa: Path, out_fa: Path, remove_ids: Set[str], dry_run: bool) -> Tuple[int, int]:
    """
    Remove matching records from a FASTA file and write the remaining records.
    """
    kept = 0
    removed = 0
    records_out = []

    for rec in SeqIO.parse(str(in_fa), "fasta"):
        if any(record_matches_id(rec, rid) for rid in remove_ids):
            removed += 1
        else:
            kept += 1
            records_out.append(rec)

    if not dry_run:
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_out, str(out_fa), "fasta")

    return kept, removed


def is_fasta(p: Path, fasta_exts: List[str]) -> bool:
    """
    Check whether a path is a FASTA file with an accepted extension.
    """
    return p.is_file() and p.suffix.lower() in set([e.lower() for e in fasta_exts])


def copy_if_missing(in_file: Path, out_dir: Path, dry_run: bool):
    """
    Copy a file to an output directory if it is not already present.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / in_file.name
    if out_file.exists():
        return
    if dry_run:
        return
    shutil.copy2(in_file, out_file)


def get_csv_list(args) -> List[Path]:
    """
    Collect frameFAIL CSV files matching the requested pattern.
    """
    csvs = list(args.fail_csv_dir.glob(args.pattern))
    if args.sort_csvs:
        csvs = sorted(csvs)
    return csvs


def main():
    """
    Run the frame-failure sequence removal workflow.
    """
    args = parse_args()

    # Determine which CSV(s) to process
    if args.single_csv:
        fail_csvs = [args.single_csv]
        full_run = False
    else:
        all_csvs = get_csv_list(args)
        if not all_csvs:
            raise SystemExit(f"No fail CSVs found with pattern '{args.pattern}' in: {args.fail_csv_dir}")

        if args.csv_index is not None:
            idx = args.csv_index - 1
            if idx < 0 or idx >= len(all_csvs):
                raise SystemExit(f"csv_index {args.csv_index} out of range (1-{len(all_csvs)})")
            fail_csvs = [all_csvs[idx]]
            full_run = False
        else:
            fail_csvs = all_csvs
            full_run = True

    if args.dry_run:
        print("DRY RUN: no outputs will be written.\n")

    # Build fail_map for the selected CSV(s)
    fail_map = {}
    for csv_path in fail_csvs:
        base_stem = csv_path.name.replace(".frameFAIL.csv", "")
        translated_ids, original_ids = load_fail_ids(csv_path, args.remove_if_passes_is)
        fail_map[base_stem] = (translated_ids, original_ids)

    print(f"Selected {len(fail_map)} fail-CSV file(s) to process.")
    if len(fail_csvs) == 1:
        print(f"CSV: {fail_csvs[0]}")

    # Process those with fail CSV
    for base_stem, (translated_remove, original_remove) in fail_map.items():
        prot_in = find_matching_fasta(base_stem, args.protein_dir, args.fasta_exts)
        cds_in = find_matching_fasta(base_stem, args.cds_dir, args.fasta_exts)

        if prot_in is None:
            print(f"[WARN] No protein FASTA found for stem: {base_stem}")
        else:
            prot_out = args.out_protein_dir / prot_in.name
            kept, removed = sift_fasta(prot_in, prot_out, translated_remove, args.dry_run)
            print(f"[PROT] {prot_in.name}: removeIDs={len(translated_remove)} removed={removed} kept={kept} -> {prot_out}")

        if cds_in is None:
            print(f"[WARN] No CDS FASTA found for stem: {base_stem}")
        else:
            cds_out = args.out_cds_dir / cds_in.name
            kept, removed = sift_fasta(cds_in, cds_out, original_remove, args.dry_run)
            print(f"[CDS ] {cds_in.name}: removeIDs={len(original_remove)} removed={removed} kept={kept} -> {cds_out}")

    # Only copy the rest in FULL RUN mode (prevents array jobs stomping each other)
    if full_run:
        for p in sorted(args.protein_dir.iterdir()):
            if not is_fasta(p, args.fasta_exts):
                continue
            if p.stem not in fail_map:
                copy_if_missing(p, args.out_protein_dir, args.dry_run)

        for p in sorted(args.cds_dir.iterdir()):
            if not is_fasta(p, args.fasta_exts):
                continue
            if p.stem not in fail_map:
                copy_if_missing(p, args.out_cds_dir, args.dry_run)

    print("\nDone.")
    print(f"SIFTED_CDS:            {args.out_cds_dir}")
    print(f"SIFTED_TRANSLATED_CDS: {args.out_protein_dir}")


if __name__ == "__main__":
    main()