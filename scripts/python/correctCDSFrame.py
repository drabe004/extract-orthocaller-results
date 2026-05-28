#!/usr/bin/env python3
import argparse
import csv
import os
import re
from typing import Dict, Tuple, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


def clean_protein(s: str) -> str:
    s = re.sub(r"\s+", "", str(s)).upper()
    while s.endswith("*"):
        s = s[:-1]
    return s


def clean_cds(s: str) -> str:
    s = re.sub(r"\s+", "", str(s)).upper()
    s = re.sub(r"[^ACGTN]", "N", s)
    return s


def species_key(seq_id: str) -> str:
    ### token up to first whitespace, then text before first underscore
    token = seq_id.split()[0]
    return token.split("_", 1)[0]


def build_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    return aligner


def align_to_gapped_strings(aligner: PairwiseAligner, s1: str, s2: str) -> Tuple[str, str]:
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


def identity_ignore_gaps(aln_a: str, aln_b: str) -> float:
    compared = 0
    matches = 0
    for aa, bb in zip(aln_a, aln_b):
        if aa == "-" or bb == "-":
            continue
        compared += 1
        if aa == bb:
            matches += 1
    return (matches / compared) if compared > 0 else 0.0


def translate_with_trim(cds: Seq, trim5: int, table: int, strand: str) -> Tuple[str, str]:
    ### Returns: (translated_protein_clean, cds_used_string)
    ### cds_used_string is after strand choice, 5' trim, and 3' trim-to-codon.
    s = clean_cds(str(cds))

    if strand == "-":
        s = str(Seq(s).reverse_complement())

    if trim5 > 0:
        if len(s) <= trim5:
            return "", ""
        s = s[trim5:]

    if len(s) < 3:
        return "", ""

    trim_len = len(s) - (len(s) % 3)
    s = s[:trim_len]

    prot = clean_protein(str(Seq(s).translate(table=table, to_stop=False)))
    return prot, s


def load_fasta_as_dict(path: str) -> Dict[str, Seq]:
    d: Dict[str, Seq] = {}
    for rec in SeqIO.parse(path, "fasta"):
        d[rec.id] = rec.seq
    return d


def build_prefix_index(ids: List[str]) -> Dict[str, List[str]]:
    ### Map prefix (before first underscore) -> list of record IDs with that prefix
    idx: Dict[str, List[str]] = {}
    for rid in ids:
        k = species_key(rid)
        idx.setdefault(k, []).append(rid)
    return idx


def resolve_protein_id_for_cds(
    cds_id: str,
    prot_dict: Dict[str, Seq],
    prot_prefix_index: Dict[str, List[str]],
) -> Tuple[Optional[str], str]:
    ### Returns (protein_id or None, match_mode)
    ### match_mode is one of: EXACT_ID, PREFIX_ID, NO_MATCH, AMBIGUOUS_PREFIX
    if cds_id in prot_dict:
        return cds_id, "EXACT_ID"

    k = species_key(cds_id)
    hits = prot_prefix_index.get(k, [])
    if len(hits) == 0:
        return None, "NO_MATCH"
    if len(hits) > 1:
        return None, "AMBIGUOUS_PREFIX"
    return hits[0], "PREFIX_ID"


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Adjust CDS (trim 0-2 nt at 5') so translation matches provided protein FASTA. "
            "Tries trim5 in {0,1,2} and optional reverse-complement. "
            "Protein matching: exact FASTA ID first, then fallback to prefix before first underscore (only if unique)."
        )
    )
    ap.add_argument("--cds_fasta", required=True, help="Input CDS fasta (one file).")
    ap.add_argument("--protein_fasta", required=True, help="Protein fasta to match (same IDs or same prefix before underscore).")
    ap.add_argument("--out_cds_fasta", required=True, help="Output corrected CDS fasta.")
    ap.add_argument("--report_csv", required=True, help="CSV report of chosen trims.")
    ap.add_argument("--genetic_code_table", type=int, default=1)
    ap.add_argument("--try_reverse_complement", action="store_true", help="Also try reverse-complement CDS.")
    ap.add_argument("--max_trim5", type=int, default=2, help="Max 5' trim to test (default 2 => tests 0,1,2).")
    ap.add_argument(
        "--require_id_match",
        action="store_true",
        help="If set, only write records where a protein can be matched; otherwise unmatched are passed through unchanged.",
    )
    args = ap.parse_args()

    cds_dict = load_fasta_as_dict(args.cds_fasta)
    prot_dict = load_fasta_as_dict(args.protein_fasta)

    prot_prefix_index = build_prefix_index(list(prot_dict.keys()))
    aligner = build_aligner()

    trim_options = list(range(0, args.max_trim5 + 1))
    strands = ["+"]
    if args.try_reverse_complement:
        strands.append("-")

    out_records = []
    rows: List[dict] = []

    n_exact_match_protein = 0
    n_changed = 0
    n_total = 0
    n_no_protein = 0
    n_ambig = 0
    n_proteins_missing_cds = 0

    cds_records = list(SeqIO.parse(args.cds_fasta, "fasta"))

    for rec in cds_records:
        n_total += 1
        cds_id = rec.id
        cds_seq = rec.seq

        prot_id, prot_match_mode = resolve_protein_id_for_cds(cds_id, prot_dict, prot_prefix_index)

        if prot_id is None:
            if prot_match_mode == "NO_MATCH":
                n_no_protein += 1
                status = "NO_PROTEIN_FOUND"
            else:
                n_ambig += 1
                status = "AMBIGUOUS_PROTEIN_PREFIX"

            if args.require_id_match:
                continue

            out_records.append(rec)
            rows.append({
                "id": cds_id,
                "status": status,
                "protein_match_mode": prot_match_mode,
                "matched_protein_id": "",
                "best_strand": "",
                "best_trim5": "",
                "best_identity": "",
                "exact_match": "",
                "orig_len_nt": len(cds_seq),
                "new_len_nt": len(cds_seq),
            })
            continue

        target_prot = clean_protein(str(prot_dict[prot_id]))

        best = {
            "identity": -1.0,
            "strand": "+",
            "trim5": 0,
            "prot": "",
            "cds_used": "",
            "exact": False,
        }

        for strand in strands:
            for trim5 in trim_options:
                cand_prot, cds_used = translate_with_trim(
                    cds_seq, trim5=trim5, table=args.genetic_code_table, strand=strand
                )
                if not cand_prot:
                    continue

                exact = (cand_prot == target_prot)
                if exact:
                    if (not best["exact"]) or (trim5 < best["trim5"]):
                        best.update({
                            "identity": 1.0,
                            "strand": strand,
                            "trim5": trim5,
                            "prot": cand_prot,
                            "cds_used": cds_used,
                            "exact": True
                        })
                    continue

                if not best["exact"]:
                    aln_a, aln_b = align_to_gapped_strings(aligner, cand_prot, target_prot)
                    ident = identity_ignore_gaps(aln_a, aln_b)
                    if ident > best["identity"]:
                        best.update({
                            "identity": ident,
                            "strand": strand,
                            "trim5": trim5,
                            "prot": cand_prot,
                            "cds_used": cds_used,
                            "exact": False
                        })

        orig_len = len(clean_cds(str(cds_seq)))

        if best["cds_used"]:
            new_seq = Seq(best["cds_used"])
            new_len = len(best["cds_used"])

            if best["exact"]:
                n_exact_match_protein += 1
            if new_len != orig_len or best["strand"] == "-" or best["trim5"] != 0:
                n_changed += 1

            out_rec = rec[:]  ### copy
            out_rec.seq = new_seq
            out_records.append(out_rec)

            rows.append({
                "id": cds_id,
                "status": "OK",
                "protein_match_mode": prot_match_mode,
                "matched_protein_id": prot_id,
                "best_strand": best["strand"],
                "best_trim5": best["trim5"],
                "best_identity": f"{best['identity']:.6f}",
                "exact_match": "YES" if best["exact"] else "NO",
                "orig_len_nt": orig_len,
                "new_len_nt": new_len,
                "target_prot_len_aa": len(target_prot),
                "translated_len_aa": len(best["prot"]),
            })
        else:
            if args.require_id_match:
                continue
            out_records.append(rec)
            rows.append({
                "id": cds_id,
                "status": "NO_VALID_TRANSLATION",
                "protein_match_mode": prot_match_mode,
                "matched_protein_id": prot_id,
                "best_strand": "",
                "best_trim5": "",
                "best_identity": "",
                "exact_match": "",
                "orig_len_nt": orig_len,
                "new_len_nt": orig_len,
            })

    for pid in prot_dict.keys():
        if pid not in cds_dict:
            n_proteins_missing_cds += 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out_cds_fasta)), exist_ok=True)
    SeqIO.write(out_records, args.out_cds_fasta, "fasta")

    os.makedirs(os.path.dirname(os.path.abspath(args.report_csv)), exist_ok=True)
    with open(args.report_csv, "w", newline="") as fh:
        if rows:
            fieldnames = list(rows[0].keys())
        else:
            fieldnames = [
                "id", "status", "protein_match_mode", "matched_protein_id", "best_strand", "best_trim5",
                "best_identity", "exact_match", "orig_len_nt", "new_len_nt"
            ]
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(
        f"DONE. total_CDS={n_total} exact_match_protein={n_exact_match_protein} changed={n_changed} "
        f"no_protein={n_no_protein} ambiguous_prefix={n_ambig} proteins_missing_CDS={n_proteins_missing_cds}"
    )


if __name__ == "__main__":
    main()
