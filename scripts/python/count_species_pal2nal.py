import os
import csv
import argparse
from Bio import SeqIO

def load_list(file):
    with open(file) as f:
        return set(line.strip() for line in f if line.strip())

def count_file(path, cavefish_set):
    species_seen = set()

    for rec in SeqIO.parse(path, "fasta"):
        sp = rec.description.split("|")[0].strip()
        species_seen.add(sp)

    total = len(species_seen)
    cave = len(species_seen & cavefish_set)
    background = total - cave

    return cave, background, total

def main():
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