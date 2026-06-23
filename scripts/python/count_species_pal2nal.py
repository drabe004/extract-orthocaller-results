```python
#!/usr/bin/env python3
"""
Count cavefish and background species represented in PAL2NAL alignments.

For each FASTA file in the input directory:
    1) Extract unique species names from FASTA headers.
    2) Count how many species belong to the supplied cavefish list.
    3) Count remaining species as background taxa.
    4) Write results to a summary CSV.

Output columns:
    file        = FASTA filename
    cavefish    = number of cavefish species present
    background  = number of non-cavefish species present
    total       = total unique species present

Typical use:
    Generate species-count summaries for filtering orthogroups prior
    to selection analyses (BUSTED, RELAX, PCOC, etc.).
"""

import os
import csv
import argparse
from Bio import SeqIO


def load_list(file):
    """
    Load a text file containing one species name per line.

    Args:
        file: Path to cavefish species list.

    Returns:
        Set of species names.
    """
    with open(file) as f:
        return set(line.strip() for line in f if line.strip())


def count_file(path, cavefish_set):
    """
    Count cavefish and background species in a FASTA file.

    Species names are extracted from the portion of each FASTA header
    before the first pipe ('|') character.

    Args:
        path: FASTA file path.
        cavefish_set: Set of cavefish species names.

    Returns:
        Tuple:
            (cavefish_count, background_count, total_species_count)
    """
    species_seen = set()

    for rec in SeqIO.parse(path, "fasta"):
        sp = rec.description.split("|")[0].strip()
        species_seen.add(sp)

    total = len(species_seen)
    cave = len(species_seen & cavefish_set)
    background = total - cave

    return cave, background, total


def main():
    """
    Parse command-line arguments and summarize species counts.

    Processes all FASTA files in the input directory matching the
    specified extensions and writes a CSV summary table.
    """
    parser = argparse.ArgumentParser(
        description="Count cavefish, background, and total species in PAL2NAL FASTA files"
    )

    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing PAL2NAL FASTA files"
    )

    parser.add_argument(
        "-c", "--cavefish_list",
        required=True,
        help="File with cavefish species names (one per line)"
    )

    parser.add_argument(
        "-o", "--output_csv",
        required=True,
        help="Output CSV file"
    )

    parser.add_argument(
        "--ext",
        nargs="+",
        default=[".fa", ".fasta", ".fas", ".fna"],
        help="FASTA file extensions (default: .fa .fasta .fas .fna)"
    )

    args = parser.parse_args()

    cavefish_set = load_list(args.cavefish_list)

    with open(args.output_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["file", "cavefish", "background", "total"])

        for fname in sorted(os.listdir(args.input_dir)):
            if not any(fname.endswith(ext) for ext in args.ext):
                continue

            path = os.path.join(args.input_dir, fname)

            cave, background, total = count_file(path, cavefish_set)

            writer.writerow([fname, cave, background, total])
            print(f"{fname}: cave={cave} background={background} total={total}")

if __name__ == "__main__":
    main()
```