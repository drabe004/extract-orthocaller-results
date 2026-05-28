#!/usr/bin/env python3
### GetOriginalProts_V1.py
### Recover unaligned/original protein sequences from per-species primary_transcripts FASTAs
### using headers from an aligned protein file. DB-aware (Ensembl / NCBI / IN HOUSE).
### Writes <input_basename>_ORIGINALSEQS.faa, per-file .nomatch.tsv, and appends a master summary TSV.

import os
import sys
import csv
import argparse
import re

### -------------------------------
### FASTA utilities
### -------------------------------

def read_fasta_iter(path):
    ### Yield (header, sequence) from a FASTA file.
    header = None
    chunks = []
    with open(path, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield (header, "".join(chunks))
                header = line.strip()[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield (header, "".join(chunks))

def index_fasta_ids(path):
    ### Index first-token IDs -> (full_header, sequence).
    idx = {}
    for h, s in read_fasta_iter(path):
        tok = h.split()[0]
        idx[tok] = (h, s)
    return idx

### -------------------------------
### Species key
### -------------------------------

def _get_first_present(row, keys):
    ### Return the first present, non-empty value among candidate keys; else raise KeyError.
    for k in keys:
        if k in row and row[k] is not None and str(row[k]).strip() != "":
            return row[k]
    raise KeyError("None of the candidate species columns present: " + ", ".join(keys))

def load_species_key_csv(csv_path):
    ### Load the species key CSV.
    ### Returns: { species_concat_no_underscores: (species_with_underscores, prot_file_name, database_type) }
    map_concat = {}
    with open(csv_path, "r", newline="") as fh:
        rdr = csv.DictReader(fh)

        species_candidates = [
            "Species Name From ... No Underscores",
            "Species",
            "Species Name From Tree"
        ]
        required = {"Prot File Name", "Database"}
        header_fields = set(rdr.fieldnames or [])
        missing_required = required - header_fields
        if missing_required:
            raise ValueError("Missing required columns in species key CSV: " + ",".join(sorted(missing_required)))
        if not any(col in header_fields for col in species_candidates):
            raise ValueError("Missing species column. Expected one of: " + ", ".join(species_candidates))

        for row in rdr:
            species_u = _get_first_present(row, species_candidates).strip()
            prot_name = row["Prot File Name"].strip()
            db_type   = row["Database"].strip()
            concat = species_u.replace("_", "")
            map_concat[concat] = (species_u, prot_name, db_type)
    return map_concat

### -------------------------------
### Header parsing
### -------------------------------

### Input headers look like:
### >TeleostClupeiformesClupeiformesClupeidaeClupeaharengus_<TOKEN>
HDR_SPLIT = re.compile(r"^>(?P<species>[A-Za-z0-9]+)_(?P<rest>\S+)")

def parse_header_line(raw_header):
    ### Return (species_concat, trailing_token) or (None, None) if not matched.
    m = HDR_SPLIT.match(">" + raw_header if not raw_header.startswith(">") else raw_header)
    if not m:
        return None, None
    return m.group("species"), m.group("rest")

### -------------------------------
### Accession candidate builders
### -------------------------------

### Ensembl-like: prefix (letters/digits) + digits, ignore trailing gene symbol
ENS_PREFIX_DIGITS = re.compile(r"^(?P<prefix>[A-Z][A-Z0-9]+?)(?P<digits>\d+).*$")
### NCBI: NP_/XP_ with/without version
NCBI_ACC = re.compile(r"(N[PX]_\d+)(\.\d+)?")
### NCBI fused: 'P0338811702...' or 'XP0338811702...' -> number + final digit is version
NCBI_FUSED = re.compile(r"^[A-Z]*P(?P<num>\d+?)(?P<ver>\d)(?!\d)")
### IN HOUSE: UN########T# -> FUN_########-T#
INHOUSE_UN = re.compile(r"UN(?P<num>\d+)T(?P<t>\d+)")

def candid_ensembl(token):
    m = ENS_PREFIX_DIGITS.match(token)
    if not m:
        return [token] if token else []
    prefix = m.group("prefix") or ""
    digits = m.group("digits") or ""
    if not digits:
        return [prefix] if prefix else []
    cands = []
    if len(digits) >= 2:
        cands.append(prefix + digits[:-1] + "." + digits[-1])
    else:
        cands.append(prefix + "0." + digits)
    cands.append(prefix + digits)
    if not prefix.startswith("ENS"):
        if len(digits) >= 2:
            cands.append("ENS" + prefix + digits[:-1] + "." + digits[-1])
        cands.append("ENS" + prefix + digits)
    seen = set(); out = []
    for c in cands:
        if c and c not in seen:
            seen.add(c); out.append(c)
    return out

def candid_ncbi(token):
    m = NCBI_ACC.search(token)
    if m:
        core = m.group(1); vers = m.group(2) or ""
        return [core + vers] if vers else [core + ".1", core]
    m2 = NCBI_FUSED.match(token)
    if m2:
        num = m2.group("num"); ver = m2.group("ver")
        return [f"NP_{num}.{ver}", f"XP_{num}.{ver}"]
    return []

def candid_inhouse(token):
    m = INHOUSE_UN.search(token)
    if not m:
        return []
    num = m.group("num"); t = m.group("t")
    return [f"FUN_{num}-T{t}"]

def derive_candidates_by_db(db_type, trailing_token):
    dbu = db_type.strip().upper()
    cands = []
    if dbu == "ENSEMBL":
        cands.extend(candid_ensembl(trailing_token))
    elif dbu == "NCBI":
        cands.extend(candid_ncbi(trailing_token))
        cands.extend(candid_ensembl(trailing_token))  ### fallback
    elif dbu == "IN HOUSE":
        cands.extend(candid_inhouse(trailing_token))
        cands.extend(candid_ensembl(trailing_token))  ### fallback
    else:
        cands.extend(candid_ncbi(trailing_token))
        cands.extend(candid_inhouse(trailing_token))
        cands.extend(candid_ensembl(trailing_token))
    seen = set(); out = []
    for c in cands:
        if c and c not in seen:
            seen.add(c); out.append(c)
    return out

### -------------------------------
### Core processing
### -------------------------------

def process_one_file(in_fasta, primary_root, key_csv, out_dir, no_match_dir, debug=False):
    species_map = load_species_key_csv(key_csv)

    base = os.path.basename(in_fasta)
    out_name = base.rsplit(".faa", 1)[0] + "_ORIGINALSEQS.faa"
    out_path = os.path.join(out_dir, out_name)

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(no_match_dir, exist_ok=True)

    idx_cache = {}
    total_in = 0
    found = 0
    misses = []
    species_files_used = set()

    with open(out_path, "w") as out_fh:
        for raw_header, _aligned_seq in read_fasta_iter(in_fasta):
            total_in += 1

            species_concat, trailing = parse_header_line(">" + raw_header if not raw_header.startswith(">") else raw_header)
            if species_concat is None:
                misses.append(("UNPARSEABLE_HEADER", raw_header, "NA"))
                if debug:
                    sys.stderr.write(f"[WARN] Unparseable header: {raw_header}\n")
                continue

            rec = species_map.get(species_concat)
            if not rec:
                fallback_key = None
                for k in species_map.keys():
                    if k.endswith(species_concat[-30:]):
                        fallback_key = k; break
                if fallback_key:
                    rec = species_map[fallback_key]
                else:
                    misses.append(("SPECIES_NOT_IN_KEY", raw_header, "NA"))
                    if debug:
                        sys.stderr.write(f"[WARN] Species not in key CSV: {species_concat}\n")
                    continue

            species_u, prot_file, db_type = rec

            ### --- minimal robust primary_transcripts path resolution ---
            base_dir = os.path.join(primary_root, "primary_transcripts")
            raw_name = prot_file.strip().rstrip("\r")
            cand_path = os.path.join(base_dir, raw_name)
            tried_paths = [cand_path]

            if not os.path.isfile(cand_path):
                if not (cand_path.endswith(".fa") or cand_path.endswith(".fasta")):
                    for ext in (".fa", ".fasta"):
                        test = cand_path + ext
                        tried_paths.append(test)
                        if os.path.isfile(test):
                            cand_path = test; break

            if not os.path.isfile(cand_path):
                name_only = os.path.basename(raw_name)
                swap = None
                if name_only.startswith("GCA_"):
                    swap = "GCF_" + name_only[4:]
                elif name_only.startswith("GCF_"):
                    swap = "GCA_" + name_only[4:]
                if swap:
                    cand2 = os.path.join(base_dir, swap)
                    tried_paths.append(cand2)
                    if not os.path.isfile(cand2) and not (cand2.endswith(".fa") or cand2.endswith(".fasta")):
                        for ext in (".fa", ".fasta"):
                            test2 = cand2 + ext
                            tried_paths.append(test2)
                            if os.path.isfile(test2):
                                cand2 = test2; break
                    if os.path.isfile(cand2):
                        cand_path = cand2

            if not os.path.isfile(cand_path):
                misses.append(("PRIMARY_FILE_MISSING", raw_header, tried_paths[0]))
                if debug:
                    sys.stderr.write(f"[ERR] Missing primary_transcripts file. Tried: {', '.join(tried_paths)}\n")
                continue

            prot_path = cand_path
            ### --- end minimal robust resolver ---

            cache_key = os.path.basename(prot_path)
            if cache_key not in idx_cache:
                if debug: sys.stderr.write(f"[INFO] Indexing {prot_path}\n")
                idx_cache[cache_key] = index_fasta_ids(prot_path)
            idx = idx_cache[cache_key]
            species_files_used.add(cache_key)

            candidates = derive_candidates_by_db(db_type, trailing)

            matched_id = None
            full_h = None
            seq = None

            ### 1) Exact first-token matches against index
            for cand in candidates:
                if cand in idx:
                    matched_id = cand
                    full_h, seq = idx[cand]
                    break

            ### 2) Fallback: partial substring scan within headers (first-token OR full header)
            if matched_id is None:
                # try each candidate as a substring
                chosen = None
                for cand in (candidates if candidates else [trailing]):
                    hits = []
                    for h, s in read_fasta_iter(prot_path):
                        first_tok = h.split()[0]
                        if (cand and (cand in first_tok or cand in h)):
                            hits.append((first_tok, h, s))
                            if len(hits) > 1:
                                break
                    if len(hits) == 1:
                        chosen = hits[0]
                        matched_id = chosen[0]
                        full_h = chosen[1]
                        seq = chosen[2]
                        break
                    elif len(hits) > 1:
                        # ambiguous; record as miss and move on to next candidate
                        if debug:
                            sys.stderr.write(f"[AMBIG] '{cand}' matched multiple headers in {prot_path}\n")
                        continue

            if matched_id is None:
                tried_cand = candidates[0] if candidates else trailing
                misses.append((tried_cand, raw_header, prot_path))
                if debug:
                    sys.stderr.write(f"[MISS] {tried_cand} in {prot_path}\n")
                continue

            ### Write out with amended header: >{species_concat}_{matched_id}
            out_h = f">{species_concat}_{matched_id}"
            out_fh.write(out_h + "\n")
            for i in range(0, len(seq), 60):
                out_fh.write(seq[i:i+60] + "\n")
            found += 1

    ### Per-file miss log
    if misses:
        nm_path = os.path.join(no_match_dir, base + ".nomatch.tsv")
        with open(nm_path, "w") as nm:
            nm.write("searched_string\tfull_input_header\tprimary_file\n")
            for s, hdr, p in misses:
                nm.write(f"{s}\t{hdr}\t{p}\n")

    ### Console summary
    sys.stdout.write("[SUMMARY] file={} total_in={} found={} misses={} species_files={}\n".format(
        base, total_in, found, len(misses), ",".join(sorted(species_files_used))
    ))

    ### Append to master summary file for easy QC aggregation
    summary_path = os.path.join(no_match_dir, "GetOriginalProts_V1.summary.tsv")
    header_needed = not os.path.exists(summary_path)
    with open(summary_path, "a") as sf:
        if header_needed:
            sf.write("file\ttotal_in\tfound\tmisses\tspecies_files\n")
        sf.write(f"{base}\t{total_in}\t{found}\t{len(misses)}\t{','.join(sorted(species_files_used))}\n")

### -------------------------------
### CLI
### -------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Recover original unaligned proteins from species primary_transcripts FASTAs (DB-aware)."
    )
    ap.add_argument("--in_fasta", required=True, help="Input aligned protein FASTA (one file)")
    ap.add_argument("--primary_dir", required=True, help="Root containing 'primary_transcripts' subdir")
    ap.add_argument(
        "--species_key_csv",
        required=True,
        help="CSV with 'Prot File Name','Species Name From ... No Underscores','Database' (accepts 'Species' or 'Species Name From Tree' as fallback)"
    )
    ap.add_argument("--out_dir", required=True, help="Output directory for recovered FASTAs")
    ap.add_argument("--no_match_dir", required=True, help="Directory for per-file .nomatch.tsv logs and master summary")
    ap.add_argument("--debug", action="store_true", help="Verbose stderr logging")
    args = ap.parse_args()

    process_one_file(
        in_fasta=args.in_fasta,
        primary_root=args.primary_dir,
        key_csv=args.species_key_csv,
        out_dir=args.out_dir,
        no_match_dir=args.no_match_dir,
        debug=args.debug
    )

if __name__ == "__main__":
    main()
