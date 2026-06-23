#!/usr/bin/env python3

"""
Translate CDS FASTA files into protein FASTA files using a specified
NCBI translation table.

The script processes all matching CDS FASTA files in an input directory,
translates nucleotide sequences to amino acid sequences, trims trailing
nucleotides when sequence lengths are not divisible by three, and writes
translated protein FASTA files to an output directory.
"""

import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    """
    Parse command-line arguments and translate CDS FASTA files
    into protein FASTA files.
    """
    parser = argparse.ArgumentParser(
        description="Translate CDS FASTA files to protein FASTA files."
    )
    parser.add_argument(
        "--indir",
        required=True,
        help="Directory containing CDS FASTA files"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to write translated protein FASTA files"
    )
    parser.add_argument(
        "--pattern",
        default="*_CDS.fasta",
        help="Glob pattern for CDS files (default: *_CDS.fasta)"
    )
    parser.add_argument(
        "--table",
        type=int,
        default=1,
        help="NCBI translation table (default: 1)"
    )

    args = parser.parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir)

    if not indir.is_dir():
        raise SystemExit(f"ERROR: {indir} is not a directory")

    outdir.mkdir(parents=True, exist_ok=True)

    cds_files = sorted(indir.glob(args.pattern))
    if not cds_files:
        raise SystemExit(f"ERROR: No files matching {args.pattern} in {indir}")

    for cds_file in cds_files:
        out_records = []

        for rec in SeqIO.parse(cds_file, "fasta"):
            seq = str(rec.seq).upper().replace(" ", "").replace("\t", "")

            trimmed_nt = 0
            if len(seq) % 3 != 0:
                trimmed_nt = len(seq) % 3
                seq = seq[:len(seq) - trimmed_nt]

            prot = Seq(seq).translate(
                table=args.table,
                to_stop=False
            )

            desc = rec.description
            if trimmed_nt > 0:
                desc = f"{desc} [trimmed_3prime_nt={trimmed_nt}]"

            out_records.append(
                SeqRecord(
                    prot,
                    id=rec.id,
                    description=desc
                )
            )

        if not out_records:
            continue

        # filename transformation
        out_name = cds_file.name
        if out_name.endswith("_CDS.fasta"):
            out_name = out_name.replace("_CDS.fasta", "_protein.faa")
        else:
            out_name = out_name + ".protein.faa"

        out_path = outdir / out_name
        SeqIO.write(out_records, out_path, "fasta")

        print(f"{cds_file.name} ? {out_path}")

    print("Done.")


if __name__ == "__main__":
    main()