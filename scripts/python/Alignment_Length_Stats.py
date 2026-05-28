import os
import sys
import csv
from statistics import mean, median

def is_fasta(fn):
    fnl = fn.lower()
    return fnl.endswith(".fa") or fnl.endswith(".fasta") or fnl.endswith(".fna")

def read_fasta(path):
    records = []
    cur_id = None
    cur_seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            records.append((cur_id, "".join(cur_seq)))
    return records

def nongap_codons(seq):
    L = len(seq)
    codons = 0
    for i in range(0, L - 2, 3):
        c = seq[i:i+3]
        if "-" in c:
            continue
        codons += 1
    return codons

def main(in_dir, out_csv, min_codons_flag):
    rows = []
    codon_lengths = []

    for fn in sorted(os.listdir(in_dir)):
        if not is_fasta(fn):
            continue

        path = os.path.join(in_dir, fn)
        recs = read_fasta(path)
        if not recs:
            continue

        nseq = len(recs)
        aln_len_nt = len(recs[0][1])

        # basic sanity: are all sequences same length?
        same_len = all(len(s) == aln_len_nt for _, s in recs)
        if not same_len:
            # still compute, but warn in the output
            pass

        aln_codons = aln_len_nt // 3 if aln_len_nt else 0

        nongap_list = [nongap_codons(seq) for _, seq in recs]
        min_ng = min(nongap_list) if nongap_list else 0
        max_ng = max(nongap_list) if nongap_list else 0
        mean_ng = mean(nongap_list) if nongap_list else 0.0

        flag_small = "YES" if aln_codons < min_codons_flag else ""

        rows.append([
            fn,
            nseq,
            aln_len_nt,
            aln_codons,
            min_ng,
            max_ng,
            round(mean_ng, 2),
            flag_small,
            "" if same_len else "WARN_not_all_same_length"
        ])

        codon_lengths.append(aln_codons)

    with open(out_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "Filename",
            "NumSeq",
            "AlnLen_nt",
            "AlnLen_codons",
            "MinNonGapCodons_perSeq",
            "MaxNonGapCodons_perSeq",
            "MeanNonGapCodons_perSeq",
            "FLAG_small_alignment",
            "Notes"
        ])
        w.writerows(rows)

    if codon_lengths:
        codon_lengths_sorted = sorted(codon_lengths)
        print("Files processed:", len(codon_lengths_sorted))
        print("Codon length min/median/max:", codon_lengths_sorted[0], int(median(codon_lengths_sorted)), codon_lengths_sorted[-1])
        small_count = sum(1 for x in codon_lengths_sorted if x < min_codons_flag)
        print("Alignments <", min_codons_flag, "codons:", small_count)
        print("Wrote:", out_csv)
    else:
        print("No FASTA files found in:", in_dir)

if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("Usage: python Alignment_Length_Stats.py <input_dir> <output_csv> [min_codons_flag]")
        sys.exit(1)

    in_dir = sys.argv[1]
    out_csv = sys.argv[2]
    min_codons_flag = int(sys.argv[3]) if len(sys.argv) == 4 else 150

    main(in_dir, out_csv, min_codons_flag)
