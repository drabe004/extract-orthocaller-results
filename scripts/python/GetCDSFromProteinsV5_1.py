#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Data import CodonTable
import hashlib, random

import argparse
import csv
import gzip
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# =========================
# Regex and helpers
# =========================
ALNUM = re.compile(r'[^0-9A-Za-z]+')
SYMBOL_FROM_FILENAME = re.compile(r'__([A-Za-z0-9\-]+)__')
ENSEMBL_GENE_ID = re.compile(r'ENS[A-Z]{2,6}G\d+(?:\.\d+)?', re.IGNORECASE)
ENSEMBL_TX_ID  = re.compile(r'ENS[A-Z]{2,6}T\d+(?:\.\d+)?', re.IGNORECASE)
ENSEMBL_VERSION = re.compile(r'\.\d+$')

# extra IDs present in CDS headers
NCBI_PROT_ID  = re.compile(r'\b(?:[A-Z]{1,3}_)?\d+\.\d+\b')   # XP_033855513.1 or 033855513.1
NCBI_BARE_ACC = re.compile(r'\b[A-Z]{2,4}\d+\.\d+\b')         # KAI2668862.1 / TNN02972.1 / etc.
FUN_ID        = re.compile(r'\bFUN[_-]?\d+-T\d+\b', re.IGNORECASE)
FUN_GARBLED   = re.compile(r'\bUN\d+T\d+\b', re.IGNORECASE)   # UN018505T1
PROTEIN_ID_FIELD = re.compile(r'\[\s*protein_id\s*=\s*([A-Z]{1,4}_[0-9]+\.[0-9]+)\s*\]', re.IGNORECASE)

CSV_SPECIES_NOUS_ALIASES = {"species name from tree no underscores"}
CSV_SPECIES_WITHUS_ALIASES = {"species name from tree"}
CSV_CDS_COL_ALIASES = {"cds"}
STD_TABLE = CodonTable.unambiguous_dna_by_id[1]  # vertebrate nuclear

# --- Ensembl-specific helpers for gene->transcript recovery via orig .pep ---
GENE_FIELD = re.compile(r'\bgene:(' + ENSEMBL_GENE_ID.pattern + r')\b', re.IGNORECASE)

def strip_ens_ver_upper(s: str) -> str:
    return strip_ens_version((s or '').upper())

def extract_ensembl_gene_from_cds_header(raw_header: str) -> Optional[str]:
    """Pull the Ensembl gene ID from a CDS header line (e.g., 'gene:ENSPNAG00000027328.1')."""
    m = GENE_FIELD.search(raw_header)
    if not m:
        # Fall back to any Ensembl gene ID present
        m2 = ENSEMBL_GENE_ID.search(raw_header)
        if not m2:
            return None
        return strip_ens_ver_upper(m2.group(0))
    return strip_ens_ver_upper(m.group(1))

def iter_fasta_filtered_by_gene(fp, wanted_gene_nov) -> Tuple[str, str]:
    """
    Yield (header, aa_seq) for records whose header contains the same Ensembl gene ID (versionless).
    We match by presence of the *versionless* gene ID substring in the header to be robust.
    """
    hdr = None
    buf = []
    for line in fp:
        if line.startswith(">"):
            if hdr is not None:
                if wanted_gene_nov in hdr.upper():
                    yield hdr, aa_clean(''.join(buf))
            hdr = line.strip()
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        if wanted_gene_nov in hdr.upper():
            yield hdr, aa_clean(''.join(buf))

def find_transcript_by_gene_and_aa_in_origpep(orig_pep_path: Path,
                                              gene_id_ver_or_not: str,
                                              target_aa: str) -> Optional[str]:
    """
    Scan the original Ensembl protein file for headers that contain the SAME gene (versionless).
    Among those, if one has AA identical to target_aa, return its Ensembl transcript ID (version STRIPPED).
    """
    if not orig_pep_path.exists() or not target_aa:
        return None
    gene_nov = strip_ens_ver_upper(gene_id_ver_or_not)
    with open_maybe_gzip(orig_pep_path) as pfh:
        candidates = []
        for pep_hdr, pep_aa in iter_fasta_filtered_by_gene(pfh, gene_nov):
            if pep_aa:
                candidates.append((pep_hdr, pep_aa))
        if not candidates:
            return None
        exact = [hdr for (hdr, aa) in candidates if aa == target_aa]
        if not exact:
            return None
        chosen_hdr = det_choice(exact, seed_key=f"{gene_nov}||{len(exact)}")
        mtx = ENSEMBL_TX_ID.search(chosen_hdr)
        if not mtx:
            return None
        return strip_ens_ver_upper(mtx.group(0))

# =========================
# mRNA fallback map (unchanged)
# =========================
MRNA_FALLBACK_MAP: Dict[str, str] = {
    "TeleostEuteleostAnabantariaSynbranchiformesMastacembelidaeMastacembelusbrichardi":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/MBRI/predict_results/Mastacembelus_brichardi.mrna-transcripts.fa",
    "TeleostOtophysaCypriniformesNemacheilidaeNemacheilustroglocateractus":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/NTROG/output_fun_NTROG1/predict_results/Nemacheilus_troglocateractus.mrna-transcripts.fa",
    "TeleostEuteleostAnabantariaSynbranchiformesSynbranchidaeOphisternoninfernale":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/OPIN/output_fun_OPIN/predict_results/Ophisternon_infernale.mrna-transcripts.fa",
    "TeleostOtophysaSiluriformesIctaluridaePrietellaphreatophila":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/Prietella/output_fun_Prietella/predict_results/Prietella_phreatophila.mrna-transcripts.fa",
    "TeleostOtophysaSiluriformesHeptapteridaeRhamdialaluchensis":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/RHLA/output_fun_RHLA/predict_results_RHLA/Rhamdia_laluchensis.mrna-transcripts.fa",
    "TeleostOtophysaSiluriformesHeptapteridaeRhamdiamacuspanensis":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/RAMA/output_fun_RAMA/predict_results/Rhamdia_macuspanensis.mrna-transcripts.fa",
    "TeleostOtophysaSiluriformesHeptapteridaeRhamdiareddelli":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/RHRE/predict_results/Rhamdia_reddelli.mrna-transcripts.fa",
    "TeleostOtophysaSiluriformesHeptapteridaeRhamdiazongolicensis":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/RHZO/output_fun_RHZO/predict_results/Rhamdia_zongolicensis.mrna-transcripts.fa",
    "TeleostOtophysaCypriniformesNemacheilidaeSchisturakaysonei":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/SKAY/predict_results/Schistura_kaysonei.mrna-transcripts.fa",
    "TeleostEuteleostOphidiariaOphidiiformesBythitidaeTyphliasanapearsei":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/TYPE/predict_results/Typhliasana_pearsei.mrna-transcripts.fa",
    "TeleostOtophysaCypriniformesNemacheilidaeSchisturaoedipus":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/SOED/predict_results/Schistura_oedipus.mrna-transcripts.fa",
    "TeleostOtophysaCypriniformesNemacheilidaeSchisturaspiesi":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/SSPI/predict_results/Schistura_spiesi.mrna-transcripts.fa",
    "TeleostOsteoglossiformesOsteoglossiformesMormyridaeStomatorhinusmicrops":
        "/projects/standard/mcgaughs/drabe004/Funnannote_Full/SMIC2/predict_results/Stomatorhinus_microps.mrna-transcripts.fa",
}

# =========================
# Utility funcs
# =========================
def translate_cds(nt: str) -> str:
    """Robust translate: allow IUPAC ambiguity, trim leftover bases, gap-tolerant.
       Ambiguous codons -> 'X'; final '*' removed if present.
    """
    if not nt:
        return ""
    s = nt.upper().replace("U", "T")
    s = re.sub(r"[\s\-\.]", "", s)  # strip whitespace and obvious gap/format chars
    aa = []
    lim = (len(s) // 3) * 3
    for i in range(0, lim, 3):
        codon = s[i:i+3]
        if re.search(r"[^ACGT]", codon):
            aa.append("X")
            continue
        try:
            aa.append(STD_TABLE.forward_table[codon])
        except KeyError:
            aa.append("*" if codon in STD_TABLE.stop_codons else "X")
    pep = "".join(aa)
    return pep[:-1] if pep.endswith("*") else pep

def aa_clean(seq: str) -> str:
    return re.sub(r'[^A-Z]', '', (seq or '').upper())

def aa_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    if L == 0:
        return 0.0
    matches = sum(1 for i in range(L) if a[i] == b[i])
    return matches / max(len(a), len(b))

def det_choice(items, seed_key: str):
    seed = int(hashlib.md5(seed_key.encode()).hexdigest(), 16) % (2**32)
    rng = random.Random(seed)
    return rng.choice(list(items))

def strip_ens_version(s: str) -> str:
    return ENSEMBL_VERSION.sub('', s)

def canonize_header(h: str) -> str:
    return re.sub(r'[\s_]+', '', h.strip().lower())

def normalize_headers(headers: List[str]) -> Dict[str, str]:
    return {canonize_header(h): h for h in headers}

def pick_header(name_map: Dict[str, str], wanted_aliases: set) -> Optional[str]:
    wanted = {canonize_header(x) for x in wanted_aliases}
    for k, v in name_map.items():
        if k in wanted:
            return v
    return None

def open_maybe_gzip(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else p.open("r", encoding="utf-8", errors="replace")

def symbol_from_filename(p: Path) -> Optional[str]:
    m = SYMBOL_FROM_FILENAME.search(p.name)
    return m.group(1) if m else None

def norm_simple(s: str) -> str:
    return ALNUM.sub('', s).upper()

# =========================
# Read CSV map (optional)
# =========================
def read_species_to_cds_map(csv_path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with csv_path.open(newline='', encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise SystemExit(f"No headers found in CSV: {csv_path}")
        name_map = normalize_headers(reader.fieldnames)
        species_nous_col = pick_header(name_map, CSV_SPECIES_NOUS_ALIASES)
        species_withus_col = pick_header(name_map, CSV_SPECIES_WITHUS_ALIASES)
        cds_col = pick_header(name_map, CSV_CDS_COL_ALIASES)
        if not cds_col:
            raise SystemExit("CSV must include a 'CDS' column.")
        if not species_nous_col and not species_withus_col:
            raise SystemExit("CSV must include 'Species Name From Tree No Underscores' or 'Species Name From Tree'.")
        for row in reader:
            cdsfile = (row.get(cds_col) or "").strip()
            if not cdsfile:
                continue
            if species_nous_col:
                key = (row.get(species_nous_col) or "").strip()
            else:
                key = (row.get(species_withus_col) or "").strip().replace("_", "")
            if key:
                mapping[key] = cdsfile
    if not mapping:
        raise SystemExit(f"No species/CDS rows parsed from CSV: {csv_path}")
    return mapping

# =========================
# Find CDS file
# =========================
def find_cds_path(cds_dir: Path, csv_value: str) -> Optional[Path]:
    candidates = []
    def add(p: Path):
        if p.exists():
            candidates.append(p)
    base = Path(csv_value).name
    add(cds_dir / base)
    add(cds_dir / (base + ".gz"))
    stem = Path(base).stem
    for ext in (".fa", ".fna"):
        add(cds_dir / (stem + ext))
        add(cds_dir / (stem + ext + ".gz"))
    return candidates[0] if candidates else None

# =========================
# Parse species + mangled
# =========================
def parse_species_prefix_and_mangled(header: str, symbol: str) -> Tuple[str, str]:
    """
    Strict rule:
      - species_key = text before the FIRST underscore `_`
      - mangled     = text AFTER that underscore and RIGHT UP TO the gene symbol
    """
    h = header[1:] if header.startswith(">") else header
    if "_" not in h:
        return h, ""
    species_key, rest = h.split("_", 1)
    if not symbol:
        return species_key, ""
    pos = rest.upper().find(symbol.upper())
    if pos == -1:
        return species_key, ""
    mangled = rest[:pos].strip()
    return species_key, mangled

# =========================
# Exact index (no fuzzy)
# =========================
class FastaExactIndex:
    __slots__ = ("path", "id_map")
    def __init__(self, path: Path):
        self.path = path
        self.id_map: Dict[str, List[str]] = {}

def add_key(idx: FastaExactIndex, key: str, raw: str):
    lst = idx.id_map.setdefault(key, [])
    if raw not in lst:
        lst.append(raw)

def index_fasta_exact(fa_path: Path) -> FastaExactIndex:
    idx = FastaExactIndex(fa_path)
    with open_maybe_gzip(fa_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            raw = line.strip()

            # Strict NCBI protein_id field
            m = PROTEIN_ID_FIELD.search(raw)
            if m:
                pid = m.group(1)
                add_key(idx, pid, raw)
                add_key(idx, norm_simple(pid), raw)

            # Ensembl IDs (version stripped)
            for m in ENSEMBL_TX_ID.finditer(raw):
                add_key(idx, strip_ens_version(m.group(0)).upper(), raw)
            for m in ENSEMBL_GENE_ID.finditer(raw):
                add_key(idx, strip_ens_version(m.group(0)).upper(), raw)

            # legacy protein_id finder (kept safe by de-dupe)
            for m in re.finditer(r'protein_id=([A-Za-z]{1,3}_\d+\.\d+)', raw):
                pid2 = m.group(1)
                add_key(idx, pid2, raw)
                add_key(idx, norm_simple(pid2), raw)

            # NCBI-like accessions
            for m in NCBI_PROT_ID.finditer(raw):
                can = m.group(0)
                add_key(idx, can, raw)
                add_key(idx, norm_simple(can), raw)

            # Bare WGS accessions
            for m in NCBI_BARE_ACC.finditer(raw):
                can = m.group(0)
                add_key(idx, can, raw)
                add_key(idx, norm_simple(can), raw)

            # FUN IDs
            for m in FUN_ID.finditer(raw):
                can = m.group(0).upper()
                add_key(idx, can, raw)
                add_key(idx, norm_simple(can), raw)

    if not idx.id_map:
        raise SystemExit(f"No headers found in FASTA: {fa_path}")
    return idx

# ---- family presence helper (limits XP/XM/NP/YP/WP guesses)
_FAM_RX = re.compile(r'^(XP|XM|NP|YP|WP)_\d+\.\d+$')

def _present_families(idx: FastaExactIndex) -> set:
    fams = set()
    for k in idx.id_map.keys():
        m = _FAM_RX.match(k)
        if m:
            fams.add(m.group(1))
    return fams

# =========================
# Reconstruction helpers (from mangled)
# =========================
def reconstruct_ncbi_prefixed(garble: str, clazz: str) -> List[str]:
    out: List[str] = []
    clazz = clazz.upper()
    last_letter = clazz[-1]
    g = norm_simple(garble)
    if g.startswith(clazz):
        tail = g[len(clazz):]
    elif g.startswith(last_letter):
        tail = g[1:]
    else:
        return out
    if not tail.isdigit():
        return out
    for ver_len in (1, 2):
        if len(tail) > ver_len:
            acc = tail[:-ver_len]
            ver = tail[-ver_len:]
            can = f"{clazz}_{acc}.{ver}"
            out.extend([can, norm_simple(can)])
    out.append(g)
    seen=set(); uniq=[]
    for x in out:
        if x not in seen:
            seen.add(x); uniq.append(x)
    return uniq

def reconstruct_ncbi_acc(garble: str) -> List[str]:
    out: List[str] = []
    G = norm_simple(garble)
    if not G or not re.match(r'^[A-Z0-9]+$', G):
        return out
    for prefix_len in (3, 2):
        if len(G) <= prefix_len + 1:
            continue
        prefix = G[:prefix_len]
        rest = G[prefix_len:]
        if not prefix.isalpha() or not rest.isdigit():
            continue
        for ver_len in (1, 2):
            if len(rest) <= ver_len:
                continue
            acc = rest[:-ver_len]
            ver = rest[-ver_len:]
            can = f"{prefix}{acc}.{ver}"
            out.extend([can, norm_simple(can)])
    out.append(G)
    seen=set(); uniq=[]
    for x in out:
        if x not in seen:
            seen.add(x); uniq.append(x)
    return uniq

def reconstruct_fun(garble: str) -> List[str]:
    t = garble.upper()
    m = re.match(r'^FUN(\d+)-T(\d+)$', t)
    if m:
        can = f"FUN{m.group(1)}-T{m.group(2)}"
        return [can, norm_simple(can)]
    m = re.match(r'^FUN(\d+)T(\d+)$', t)
    if m:
        can = f"FUN{m.group(1)}-T{m.group(2)}"
        return [can, norm_simple(can)]
    m = re.match(r'^UN(\d+)T(\d+)$', t)
    if m:
        can = f"FUN{m.group(1)}-T{m.group(2)}"
        return [can, norm_simple(can), t, norm_simple(t)]
    return [norm_simple(t)]

def rebuilt_candidates_from_mangled(mangled: str, scheme: Dict[str, Any]) -> List[str]:
    cands: List[str] = []
    stype = (scheme.get("type") or "").lower()
    clazz = (scheme.get("class") or "").upper()
    M = (mangled or "").upper()
    if not M:
        return []
    if stype == "ncbi":
        if clazz in {"XP", "XM", "NP", "YP", "WP"}:
            toks = re.findall(r'(?:[A-Z]{1,3}_?\d{6,}(?:\.\d{1,2})?|\d{6,}(?:\.\d{1,2})?)', M)
            if not toks:
                toks = [M]
            single_to_family = {"P": "XP", "N": "NP", "Y": "YP", "W": "WP"}
            family_priority  = ("XP", "NP")
            for tok in toks:
                T = tok.upper()
                m = re.match(r'^([A-Z]{1,3})?_?(\d{6,})(?:\.(\d{1,2}))?$', T)
                if not m:
                    cands.extend(reconstruct_ncbi_prefixed(T, clazz))
                    continue
                fam_hint = m.group(1)
                acc_part = m.group(2)
                ver_dot  = m.group(3)
                if ver_dot:
                    acc = acc_part; ver = ver_dot
                else:
                    if len(acc_part) < 2:
                        continue
                    acc = acc_part[:-1]; ver = acc_part[-1]
                bare = f"{acc}.{ver}"
                cands.append(bare)
                fams: List[str] = []
                if fam_hint:
                    fams.append(single_to_family.get(fam_hint, fam_hint))
                else:
                    fams.extend(family_priority)
                for fam in fams:
                    cands.append(f"{fam}_{acc}.{ver}")
                cands.append(T)
        else:
            runs = re.findall(r'[A-Z0-9]{5,}', M)
            if not runs:
                runs = [M]
            for r in runs:
                cands.extend(reconstruct_ncbi_acc(r))
    elif stype == "fun":
        runs = re.findall(r'(?:FUN[_-]?\d+-?T\d+|UN\d+T\d+)', M)
        if not runs:
            runs = [M]
        for r in runs:
            cands.extend(reconstruct_fun(r))
    elif stype == "ens_only":
        m_gene = re.search(r'([A-Z]{2,12}G)(\d+)(?:\.(\d+))?', M, flags=re.IGNORECASE)
        if m_gene:
            head = m_gene.group(1).upper()
            digits = m_gene.group(2)
            ver_dot = m_gene.group(3)
            tok_nodot = f"{head}{digits}"
            cands.append(tok_nodot)
            if not ver_dot and len(digits) >= 2:
                base1, ver1 = digits[:-1], digits[-1]
                tok_v1 = f"{head}{base1}.{ver1}"
                cands.append(tok_v1)
                if len(digits) >= 3:
                    base2, ver2 = digits[:-2], digits[-2:]
                    tok_v2 = f"{head}{base2}.{ver2}"
                    cands.append(tok_v2)
            if head.startswith("NS"):
                ens_head = "ENS" + head[2:]
                cands.append(f"{ens_head}{digits}")
                if ver_dot:
                    cands.append(f"{ens_head}{digits}.{ver_dot}")
                else:
                    if len(digits) >= 2:
                        cands.append(f"{ens_head}{digits[:-1]}.{digits[-1]}")
                        if len(digits) >= 3:
                            cands.append(f"{ens_head}{digits[:-2]}.{digits[-2:]}")
    seen = set(); uniq: List[str] = []
    for k in cands:
        if k not in seen:
            seen.add(k); uniq.append(k)
    return uniq

# ---- conservative rescues ----
def _first_letter_rescue(idx: FastaExactIndex, candidate: str) -> List[str]:
    if '.' not in candidate or len(candidate) < 6:
        return []
    protein_keys = [k for k in idx.id_map.keys() if re.match(r'^[A-Z]{3}\d+\.\d+$', k)]
    matches = [k for k in protein_keys if k[1:] == candidate]
    if len(matches) == 1:
        return list(idx.id_map[matches[0]])
    return []

def _unique_substring_fallback(idx: FastaExactIndex, candidates: List[str]) -> Tuple[List[str], Optional[str]]:
    if not candidates:
        return [], None
    keys = list(idx.id_map.keys())
    for c in candidates:
        if not c or '.' not in c or len(c) < 6:
            continue
        matches = [k for k in keys if c in k]
        if len(matches) == 1:
            k = matches[0]
            return list(idx.id_map[k]), k
    return [], None

def try_exact_match(idx: FastaExactIndex,
                    prot_header: str,
                    mangled: str,
                    scheme: Dict[str, Any],
                    dbg_emitter=None) -> Tuple[List[str], List[str]]:
    tried: List[str] = []
    can_keys: List[str] = []
    norm_keys: List[str] = []
    for key in rebuilt_candidates_from_mangled(mangled, scheme):
        if not key:
            continue
        can_keys.append(key)
        if key.upper().startswith("ENS"):
            stripped = strip_ens_version(key.upper())
            can_keys.append(stripped)
            norm_keys.append(stripped)
        else:
            norm_keys.append(norm_simple(key))
    def _uniq(xs):
        seen=set(); out=[]
        for x in xs:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    can_keys  = _uniq(can_keys)
    norm_keys = _uniq(norm_keys)

    def _reject_trailing_ref(k: str) -> bool:
        return not k.upper().startswith("ENSDARG")  # avoid zebrafish IDs in other spp
    can_keys  = [k for k in can_keys  if _reject_trailing_ref(k)]
    norm_keys = [k for k in norm_keys if _reject_trailing_ref(k)]

    fams_here = _present_families(idx)
    def _family_ok(key: str) -> bool:
        m = re.match(r'^(XP|XM|NP|YP|WP)_', key)
        return (m is None) or (m.group(1) in fams_here)
    can_keys  = [k for k in can_keys  if _family_ok(k)]
    norm_keys = [k for k in norm_keys if _family_ok(k)]

    for k in can_keys:
        tried.append(k)
        if k in idx.id_map:
            return list(idx.id_map[k]), tried
    for k in norm_keys:
        tried.append(k)
        if k in idx.id_map:
            return list(idx.id_map[k]), tried
    for c in can_keys:
        hits = _first_letter_rescue(idx, c)
        if hits:
            tried.append(f"[1stLETTER]{c}")
            return hits, tried
    hits, matched_key = _unique_substring_fallback(idx, can_keys)
    if hits:
        tried.append(f"[SUBSTR]{matched_key}")
        return hits, tried
    return [], tried

# =========================
# IO helpers
# =========================
def extract_sequence_at_header(fa_path: Path, target_header: str) -> str:
    seq_lines = []
    capture = False
    with open_maybe_gzip(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if capture:
                    break
                capture = (line.strip() == target_header)
                continue
            if capture:
                seq_lines.append(line.strip())
    return "".join(seq_lines)

def seq_len_for_header(fa_path: Path, target_header: str) -> int:
    return len(extract_sequence_at_header(fa_path, target_header))

def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def out_name_for(protein_faa: Path, out_dir: Optional[Path], suffix: str = "_CDS.fasta") -> Path:
    name = protein_faa.name + suffix
    return (out_dir / name) if out_dir else protein_faa.with_name(name)

def ambcalls_name_for(protein_faa: Path, ambcalls_dir: Optional[Path], out_dir: Optional[Path], suffix: str = ".ambcalls.tsv") -> Path:
    base_dir = ambcalls_dir or out_dir or protein_faa.parent
    return base_dir / (protein_faa.name + suffix)

def iter_fasta(fp):
    hdr = None
    buf = []
    for line in fp:
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, ''.join(buf)
            hdr = line.strip()
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        yield hdr, ''.join(buf)

# =========================
# Ensembl orig-protein indexing (NEW)
# =========================
def extract_ids_from_header_text(raw_header: str) -> Tuple[List[str], List[str]]:
    genes = [strip_ens_version(m.group(0)).upper() for m in ENSEMBL_GENE_ID.finditer(raw_header)]
    txs   = [strip_ens_version(m.group(0)).upper() for m in ENSEMBL_TX_ID.finditer(raw_header)]
    # dedup preserving order
    seen=set(); genes_u=[]; txs_u=[]
    for g in genes:
        if g not in seen:
            seen.add(g); genes_u.append(g)
    seen.clear()
    for t in txs:
        if t not in seen:
            seen.add(t); txs_u.append(t)
    return genes_u, txs_u

def _species_key_variants_from_filename(fname_stem: str) -> Tuple[str, str]:
    """
    For a file like 'Amphilophus_citrinellus.Midas_v5.pep.all.fa',
    this returns ('Amphilophus_citrinellus', 'Amphilophuscitrinellus').
    """
    sp = fname_stem.split('.', 1)[0]
    return sp, sp.replace('_','')

class OrigProtIndex:
    __slots__ = ("path", "by_aa_md5", "by_gene", "by_tx", "header_to_seq")
    def __init__(self, path: Path):
        self.path = path
        self.by_aa_md5: Dict[str, List[str]] = {}
        self.by_gene: Dict[str, List[str]] = {}
        self.by_tx:   Dict[str, List[str]] = {}
        self.header_to_seq: Dict[str, str] = {}

def _add_map_list(d: Dict[str, List[str]], k: str, h: str):
    if not k:
        return
    lst = d.setdefault(k, [])
    if h not in lst:
        lst.append(h)

def index_original_proteins(prot_path: Path) -> OrigProtIndex:
    opi = OrigProtIndex(prot_path)
    with open_maybe_gzip(prot_path) as fh:
        for hdr, seq in iter_fasta(fh):
            if not hdr.startswith(">"):
                continue
            aa = aa_clean(seq)
            if not aa:
                continue
            md5 = hashlib.md5(aa.encode()).hexdigest()
            opi.header_to_seq[hdr] = aa
            _add_map_list(opi.by_aa_md5, md5, hdr)
            genes, txs = extract_ids_from_header_text(hdr)
            for g in genes:
                _add_map_list(opi.by_gene, g, hdr)
            for t in txs:
                _add_map_list(opi.by_tx, t, hdr)
    if not opi.header_to_seq:
        raise SystemExit(f"No protein entries indexed in: {prot_path}")
    return opi

# Global caches (NEW)
ORIG_PROT_CACHE: Dict[str, OrigProtIndex] = {}   # key = absolute prot file path
SPECIES_TO_ORIGPROT: Dict[str, Path] = {}        # maps both 'Genus_species' and 'Genusspecies' -> Path

def build_species_to_origprot(orig_dir: Path):
    """
    Scan --orig_proteins_dir for files like *.pep.all.fa[.gz]
    Index under both 'Genus_species' and 'Genusspecies' keys.
    """
    patterns = ["*.pep.all.fa", "*.pep.all.fa.gz", "*.fa", "*.fa.gz"]
    for pat in patterns:
        for p in orig_dir.glob(pat):
            stem_first = p.name.split(".fa")[0]  # tolerate .fa or .fa.gz
            species_with_us, species_no_us = _species_key_variants_from_filename(stem_first)
            SPECIES_TO_ORIGPROT.setdefault(species_with_us, p)
            SPECIES_TO_ORIGPROT.setdefault(species_no_us, p)

def filter_ensembl_gene_candidates(keys: List[str]) -> List[str]:
    out = []
    seen = set()
    for k in keys:
        K = strip_ens_version(k.upper())
        if ENSEMBL_GENE_ID.fullmatch(K) and K not in seen:
            seen.add(K)
            out.append(K)
    return out

def resolve_with_original_proteins(
    species_key: str,
    prot_seq: str,
    ens_gene_candidates: List[str],
    cds_idx: FastaExactIndex,
    cds_hits: List[str],
) -> Optional[str]:
    """
    Disambiguate Ensembl multi-hits via original protein FASTA for this species:
      - restrict to original headers carrying one of ens_gene_candidates
      - require identical AA to prot_seq
      - take transcript ID(s) from those headers
      - intersect with current CDS hits via cds_idx.id_map
    """
    if not prot_seq:
        return None
    prot_path = SPECIES_TO_ORIGPROT.get(species_key) or SPECIES_TO_ORIGPROT.get(species_key.replace("_",""))
    if not prot_path or not prot_path.exists():
        return None

    cache_key = str(prot_path.resolve())
    if cache_key not in ORIG_PROT_CACHE:
        ORIG_PROT_CACHE[cache_key] = index_original_proteins(prot_path)
    opi = ORIG_PROT_CACHE[cache_key]

    md5 = hashlib.md5(aa_clean(prot_seq).encode()).hexdigest()

    # 1) restrict to headers with candidate genes
    gene_keys = {strip_ens_version(g.upper()) for g in ens_gene_candidates}
    candidate_headers = set()
    for g in gene_keys:
        for h in opi.by_gene.get(g, []):
            candidate_headers.add(h)

    # Fallback: if no gene-match, at least require identical AA across all originals
    if not candidate_headers:
        candidate_headers = set(opi.by_aa_md5.get(md5, []))
    else:
        candidate_headers = {h for h in candidate_headers if hashlib.md5(opi.header_to_seq.get(h,"").encode()).hexdigest() == md5}

    if not candidate_headers:
        return None

    # 2) collect transcript IDs from those headers
    tx_ids = []
    seen_t = set()
    for h in sorted(candidate_headers):
        _g, txs = extract_ids_from_header_text(h)
        for t in txs:
            T = strip_ens_ver_upper(t)
            if T.startswith("ENS") and T not in seen_t:
                seen_t.add(T)
                tx_ids.append(T)
    if not tx_ids:
        return None

    # 3) map transcript IDs into current cds_hits via the CDS index
    tx_hit_headers = []
    for t in tx_ids:
        for raw_h in cds_idx.id_map.get(t, []):
            if raw_h in cds_hits:
                tx_hit_headers.append(raw_h)

    tx_hit_headers = sorted(set(tx_hit_headers))
    if not tx_hit_headers:
        return None
    return tx_hit_headers[0]  # deterministic

# =========================
# Core processing
# =========================
def resolve_by_aa_len_then_detpick(
    target_aa: str,
    candidates,           # list of dicts: {'header','cds_seq','protein_id','gene','source_file'}
    og_key: str,
    prot_header: str,
):
    for c in candidates:
        if 'aa_seq' not in c or c['aa_seq'] is None:
            nt = c.get('cds_seq', '')
            c['aa_seq'] = translate_cds(nt) if nt else ''
    aa_exact = [c for c in candidates if target_aa and c['aa_seq'] == target_aa]
    if aa_exact:
        pool = aa_exact; reason_base = "AA_EXACT"
    else:
        pool = candidates; reason_base = None
    max_len = max((len(c['aa_seq']) for c in pool), default=-1)
    len_pool = [c for c in pool if len(c['aa_seq']) == max_len]
    if len(len_pool) == 1:
        reason = reason_base or "AA_LEN_MAX"
        return len_pool[0], reason, 0, []
    chosen = det_choice(len_pool, seed_key=f"{og_key}||{prot_header}")
    others = [c for c in len_pool if c is not chosen]
    return chosen, (reason_base or "DET_PICK"), len(len_pool)-1, others

# --- Ambiguity logging helpers (NEW) ---
def _md5(s: str) -> str:
    return hashlib.md5((s or "").encode()).hexdigest()

def _join(xs):
    return ";".join(xs) if xs else ""

def _lens(xs):
    return ";".join(str(x) for x in xs) if xs else ""

def log_ambcall(
    recs_list: List[str],
    species_key: str,
    file_name: str,
    protein_header: str,
    chosen_header: str,
    chosen_len: int,
    others: List[str],
    other_lens: List[int],
    reason: str,
    tie_path: List[str],
    query_aa_seq: str,
    orig_aa_seq: Optional[str],
    all_headers: List[str],
    all_aa_lens: List[int],
):
    """
    Appends a single, extended ambcalls TSV row.
    Columns (added after your existing ones):
      tie_path      -> '>'-separated steps showing how the tie was resolved
      query_aa      -> AA sequence of the query protein (as seen in the proteins_fasta)
      orig_aa       -> AA sequence from original Ensembl protein (when used; else empty)
      all_headers   -> ALL ambiguous CDS headers considered (chosen + others)
      all_aa_lens   -> AA lengths for ALL headers in all_headers (same order)
    """
    recs_list.append("\t".join([
        species_key, file_name, protein_header,
        chosen_header, str(chosen_len), str(len(others)),
        _join(others), _lens(other_lens),
        reason,
        ">".join(tie_path) if tie_path else "",
        query_aa_seq or "",
        (orig_aa_seq or ""),
        _join(all_headers),
        _lens(all_aa_lens),
    ]))

def process_one_file(
    proteins_path: Path,
    cds_dir: Path,
    out_dir: Optional[Path],
    ambcalls_dir: Optional[Path] = None,
    no_match_log_dir: Optional[Path] = None,
    primary_tie_break: bool = False,
    scheme_lookup: Optional[Dict[str, Dict[str, Any]]] = None,
    species_to_cds: Optional[Dict[str, str]] = None,
    debug: bool = False,
    orig_proteins_dir: Optional[Path] = None,   # <-- NEW
) -> Tuple[int,int,int,Optional[Path],Optional[Path],Optional[Path]]:

    scheme_lookup = scheme_lookup or {}
    species_to_cds = species_to_cds or {}

    symbol = symbol_from_filename(proteins_path)
    out_path = out_name_for(proteins_path, out_dir)

    if not symbol:
        nomatch_path = (no_match_log_dir or out_dir or proteins_path.parent) / (out_path.name + ".nomatch.txt")
        nomatch_path.parent.mkdir(parents=True, exist_ok=True)
        with open(nomatch_path, "w", encoding="utf-8") as nm:
            nm.write("TAG\tFILE\tPROT_HEADER\tEXTRA\n")
            nm.write(f"missing_symbol_in_filename\tNA\t{proteins_path.name}\tneed pattern like '__ACKR2__'\n")
        return 0, 0, 0, out_path, None, nomatch_path

    ambcalls_path = ambcalls_name_for(proteins_path, ambcalls_dir, out_dir)
    nomatch_path = (no_match_log_dir or out_dir or proteins_path.parent) / (out_path.name + ".nomatch.txt")

    index_cache: Dict[str, FastaExactIndex] = {}
    total = 0
    unamb_matched = 0
    resolved_ambig = 0
    nomatch_records: List[str] = []
    ambcalls_records: List[str] = []

    with open(proteins_path, "r", encoding="utf-8", errors="replace") as pfh, \
         open(out_path, "w", encoding="utf-8") as ofh:

        for prot_header, prot_seq_raw in iter_fasta(pfh):
            total += 1
            prot_seq = aa_clean(prot_seq_raw)

            species_key, mangled = parse_species_prefix_and_mangled(prot_header, symbol)

            cds_file_in_map = species_to_cds.get(species_key)
            if not cds_file_in_map:
                nomatch_records.append(f"no matches\tNA\t{prot_header}\tNo CDS file for species_key={species_key}")
                continue

            cds_path = find_cds_path(cds_dir, cds_file_in_map)
            if not cds_path:
                nomatch_records.append(f"no matches\t{cds_file_in_map}\t{prot_header}\tCDS file not found")
                continue

            key = str(cds_path.resolve())
            if key not in index_cache:
                index_cache[key] = index_fasta_exact(cds_path)
            idx = index_cache[key]

            scheme = scheme_lookup.get(species_key, {"type": "ens_only"})
            hits, tried_keys = try_exact_match(idx, prot_header, mangled, scheme)

            if debug:
                print(f"[DBG] species_key={species_key}  header={prot_header}")
                print(f"[DBG] cds_file={cds_file_in_map} -> {cds_path.name}")
                print(f"[DBG] scheme={scheme}")
                show = ", ".join(tried_keys[:8])
                print(f"[DBG] tried_keys[:8]={show if show else '[]'}")
                if not hits:
                    print(f"[DBG] no exact hits in CDS")

            # strict mRNA fallback for listed DIY genomes (exact logic)
            if not hits and species_key in MRNA_FALLBACK_MAP:
                mrna_path = Path(MRNA_FALLBACK_MAP[species_key])
                if mrna_path.exists():
                    mkey = str(mrna_path.resolve())
                    if mkey not in index_cache:
                        index_cache[mkey] = index_fasta_exact(mrna_path)
                    midx = index_cache[mkey]
                    hits, tried_keys2 = try_exact_match(midx, prot_header, mangled, scheme)
                    tried_keys.extend(tried_keys2)

                    if hits:
                        candidates = []
                        for hdr in hits:
                            nt_seq = extract_sequence_at_header(mrna_path, hdr)
                            candidates.append({'header': hdr, 'cds_seq': nt_seq})

                        chosen, reason, n_tied, others = resolve_by_aa_len_then_detpick(
                            target_aa=prot_seq,
                            candidates=candidates,
                            og_key=species_key,
                            prot_header=prot_header,
                        )

                        # extended ambcalls
                        tie_path = ["mrna_fallback", "aa_exact_or_len_detpick"]
                        all_headers = [c['header'] for c in ([chosen] + others)]
                        all_aa_lens = [len((c.get('aa_seq') or "")) for c in ([chosen] + others)]
                        other_hdrs = [o['header'] for o in others]
                        other_lens = [len((o.get('cds_seq') or "")) for o in others]

                        log_ambcall(
                            ambcalls_records,
                            species_key,
                            mrna_path.name,
                            prot_header,
                            chosen['header'],
                            len(chosen.get('cds_seq') or ""),
                            other_hdrs,
                            other_lens,
                            "multiple hits",
                            tie_path,
                            prot_seq,
                            None,
                            all_headers,
                            all_aa_lens,
                        )

                        raw_header = chosen['header']
                        seq = chosen['cds_seq'] or extract_sequence_at_header(mrna_path, raw_header)
                        if not seq:
                            nomatch_records.append(
                                f"no matches\t{mrna_path.name}\t{prot_header}\tMatched mRNA header but no sequence: {raw_header}"
                            )
                            continue

                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue  # done with this protein entry

            if not hits:
                nomatch_records.append(f"no matches\t{cds_path.name}\t{prot_header}\tCANDS={';'.join(tried_keys) if tried_keys else 'NONE'}")
                continue

            # --- Ensembl ambiguity resolution via original .pep: gene -> AA-identical -> transcript ID ---
            if len(hits) > 1 and (scheme.get("type") or "").lower() == "ens_only" and prot_seq:
                # Extract (versionless) Ensembl gene IDs from the ambiguous CDS headers
                gene_ids = []
                for h in hits:
                    gid = extract_ensembl_gene_from_cds_header(h)
                    if gid:
                        gene_ids.append(gid)
                gene_ids = list({g for g in gene_ids})  # unique

                chosen_from_orig = None

                # Preferred path: use cached species?orig-protein index
                if gene_ids and SPECIES_TO_ORIGPROT:
                    chosen_hdr = resolve_with_original_proteins(
                        species_key=species_key,
                        prot_seq=prot_seq,
                        ens_gene_candidates=gene_ids,
                        cds_idx=idx,
                        cds_hits=hits,
                    )
                    if chosen_hdr:
                        chosen_from_orig = chosen_hdr

                # Fallback: derive .pep filename from the CDS filename if an orig_proteins_dir was provided
                if not chosen_from_orig and gene_ids and orig_proteins_dir:
                    cds_base = Path(cds_path).name
                    # *.cds.all.fa(.gz) -> *.pep.all.fa(.gz)
                    pep_guess = re.sub(r'\.cds\.all\.fa(\.gz)?$', r'.pep.all.fa\1', cds_base)
                    guess_path = orig_proteins_dir / pep_guess
                    if guess_path.exists():
                        for gid in gene_ids:
                            tx_id = find_transcript_by_gene_and_aa_in_origpep(guess_path, gid, prot_seq)
                            if tx_id:
                                tx_nov = strip_ens_ver_upper(tx_id)
                                preferred = [h for h in hits if tx_nov in h.upper()]
                                if preferred:
                                    chosen_from_orig = sorted(preferred)[0]
                                    break

                if chosen_from_orig:
                    raw_header = chosen_from_orig
                    seq = extract_sequence_at_header(cds_path, raw_header)
                    if seq:
                        if debug:
                            print(f"[DBG] resolved via orig-prot to transcript in header: {raw_header}")

                        # extended ambcalls for orig-prot path
                        tie_path = ["orig_prot", "gene_filtered", "aa_exact", "tx_match_in_cds"]
                        orig_aa_seq = prot_seq  # identical by construction here

                        all_headers = sorted(set(hits))
                        all_aa_lens = []
                        for h in all_headers:
                            aa_tmp = translate_cds(extract_sequence_at_header(cds_path, h))
                            all_aa_lens.append(len(aa_tmp or ""))

                        others = sorted(h for h in hits if h != raw_header)
                        other_lens = [len(extract_sequence_at_header(cds_path, h)) for h in others]

                        log_ambcall(
                            ambcalls_records,
                            species_key,
                            cds_path.name,
                            prot_header,
                            raw_header,
                            len(seq),
                            others,
                            other_lens,
                            "aa_exact_from_origprot_gene_filtered",
                            tie_path,
                            prot_seq,
                            orig_aa_seq,
                            all_headers,
                            all_aa_lens,
                        )

                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue

            # --- AA-based disambiguation before generic tie-breaks ---
            if len(hits) > 1 and prot_seq:
                translations = []
                for h in hits:
                    aa = translate_cds(extract_sequence_at_header(cds_path, h))
                    translations.append((h, aa))

                exact = [h for (h, aa) in translations if aa and aa == prot_seq]
                if len(exact) == 1:
                    raw_header = exact[0]
                    seq = extract_sequence_at_header(cds_path, raw_header)
                    if not seq:
                        nomatch_records.append(f"no matches\t{cds_path.name}\t{prot_header}\tMatched by AA but no sequence: {raw_header}")
                    else:
                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue
                elif len(exact) > 1:
                    hits = exact

                scored = []
                for h, aa in translations:
                    if not aa:
                        continue
                    scored.append((aa_identity(prot_seq, aa), h))
                if scored:
                    scored.sort(key=lambda x: (-x[0], x[1]))
                    best_id, best_h = scored[0]
                    if len(scored) == 1 or best_id > scored[1][0]:
                        raw_header = best_h
                        seq = extract_sequence_at_header(cds_path, raw_header)
                        if not seq:
                            nomatch_records.append(f"no matches\t{cds_path.name}\t{prot_header}\tMatched by AA identity but no sequence: {raw_header}")
                        else:
                            out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                            ofh.write(out_header + "\n")
                            ofh.write(wrap(seq) + "\n")
                            unamb_matched += 1
                            resolved_ambig += 1
                            continue

            # --- exact-AA/closest-length/NT-identity resolution ---
            if len(hits) > 1 and prot_seq:
                nt_by_hdr = {}
                aa_by_hdr = {}
                aa_len_by_hdr = {}
                for h in hits:
                    nt = extract_sequence_at_header(cds_path, h)
                    nt_by_hdr[h] = nt
                    aa = translate_cds(nt)
                    aa_by_hdr[h] = aa
                    aa_len_by_hdr[h] = len(aa)

                prot_len = len(prot_seq)
                exact_aa_hits = [h for h in hits if aa_by_hdr.get(h) == prot_seq]

                if len(exact_aa_hits) == 1:
                    raw_header = exact_aa_hits[0]
                    seq = nt_by_hdr[raw_header]
                    if not seq:
                        nomatch_records.append(f"no matches\t{cds_path.name}\t{prot_header}\tMatched by exact AA but no sequence: {raw_header}")
                    else:
                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue

                elif len(exact_aa_hits) > 1:
                    deltas = [(abs(aa_len_by_hdr[h] - prot_len), h) for h in exact_aa_hits]
                    min_delta = min(deltas, key=lambda x: x[0])[0]
                    closest = [h for d, h in deltas if d == min_delta]

                    if len(closest) == 1:
                        raw_header = closest[0]
                        seq = nt_by_hdr[raw_header]
                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue

                    nts = {nt_by_hdr[h] for h in closest}
                    if len(nts) == 1:
                        chosen = sorted(closest)[0]
                        chosen_len = len(next(iter(nts)))
                        others = sorted(h for h in closest if h != chosen)
                        other_lens = [len(nt_by_hdr[h]) for h in others]

                        # extended ambcalls
                        tie_path = ["cds", "aa_exact_multi", "closest_len", "nt_identical", "det_pick"]
                        all_headers = sorted(set(closest))
                        all_aa_lens = []
                        for h2 in all_headers:
                            aa_tmp = translate_cds(nt_by_hdr[h2])
                            all_aa_lens.append(len(aa_tmp or ""))

                        log_ambcall(
                            ambcalls_records,
                            species_key,
                            cds_path.name,
                            prot_header,
                            chosen,
                            chosen_len,
                            others,
                            other_lens,
                            "multiple hits",
                            tie_path,
                            prot_seq,
                            None,
                            all_headers,
                            all_aa_lens,
                        )

                        raw_header = chosen
                        seq = nt_by_hdr[chosen]
                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue
                    else:
                        chosen = sorted(closest)[0]
                        chosen_len = len(nt_by_hdr[chosen])
                        others = sorted(h for h in closest if h != chosen)
                        other_lens = [len(nt_by_hdr[h]) for h in others]

                        # extended ambcalls
                        tie_path = ["cds", "aa_exact_multi", "closest_len", "det_pick"]
                        all_headers = sorted(set(closest))
                        all_aa_lens = []
                        for h2 in all_headers:
                            aa_tmp = translate_cds(nt_by_hdr[h2])
                            all_aa_lens.append(len(aa_tmp or ""))

                        log_ambcall(
                            ambcalls_records,
                            species_key,
                            cds_path.name,
                            prot_header,
                            chosen,
                            chosen_len,
                            others,
                            other_lens,
                            "multiple hits",
                            tie_path,
                            prot_seq,
                            None,
                            all_headers,
                            all_aa_lens,
                        )

                        raw_header = chosen
                        seq = nt_by_hdr[chosen]
                        out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                        ofh.write(out_header + "\n")
                        ofh.write(wrap(seq) + "\n")
                        unamb_matched += 1
                        resolved_ambig += 1
                        continue

            # --- collapse identical sequences, then pick deterministically if still tied ---
            if len(hits) > 1:
                nt_by_hdr = {h: extract_sequence_at_header(cds_path, h) for h in hits}
                uniq_nt_groups: Dict[str, List[str]] = {}
                for h, s in nt_by_hdr.items():
                    uniq_nt_groups.setdefault(s, []).append(h)
                if len(uniq_nt_groups) == 1:
                    chosen = sorted(hits)[0]
                    chosen_len = len(next(iter(nt_by_hdr.values())))
                    others = sorted(h for h in hits if h != chosen)
                    other_lens = [len(nt_by_hdr[h]) for h in others]

                    # extended ambcalls
                    tie_path = ["cds", "nt_identical_collapse", "det_pick"]
                    all_headers = sorted(set(hits))
                    all_aa_lens = []
                    for h2 in all_headers:
                        aa_tmp = translate_cds(nt_by_hdr[h2])
                        all_aa_lens.append(len(aa_tmp or ""))

                    log_ambcall(
                        ambcalls_records,
                        species_key,
                        cds_path.name,
                        prot_header,
                        chosen,
                        chosen_len,
                        others,
                        other_lens,
                        "multiple hits",
                        tie_path,
                        prot_seq,
                        None,
                        all_headers,
                        all_aa_lens,
                    )

                    raw_header = chosen
                    seq = nt_by_hdr[chosen]
                    out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                    ofh.write(out_header + "\n")
                    ofh.write(wrap(seq) + "\n")
                    unamb_matched += 1
                    resolved_ambig += 1
                    continue

            # Final tie-break
            if len(hits) > 1:
                candidates = []
                for hdr in hits:
                    nt_seq = extract_sequence_at_header(cds_path, hdr)
                    candidates.append({'header': hdr, 'cds_seq': nt_seq})
                chosen, reason, n_tied, others = resolve_by_aa_len_then_detpick(
                    target_aa=prot_seq,
                    candidates=candidates,
                    og_key=species_key,
                    prot_header=prot_header,
                )
                chosen_len_nt = len(chosen['cds_seq'] or "")
                other_hdrs = [o['header'] for o in others]
                other_lens = [len(o['cds_seq'] or "") for o in others]

                # extended ambcalls
                tie_path = ["cds", "aa_len_max" if reason != "AA_EXACT" else "aa_exact", "det_pick"]
                all_headers = [c['header'] for c in candidates]
                all_aa_lens = []
                for c in candidates:
                    aa_tmp = c.get('aa_seq')
                    if aa_tmp is None:
                        aa_tmp = translate_cds(c.get('cds_seq') or "")
                    all_aa_lens.append(len(aa_tmp or ""))

                log_ambcall(
                    ambcalls_records,
                    species_key,
                    cds_path.name,
                    prot_header,
                    chosen['header'],
                    chosen_len_nt,
                    other_hdrs,
                    other_lens,
                    "multiple hits",
                    tie_path,
                    prot_seq,
                    None,
                    all_headers,
                    all_aa_lens,
                )

                raw_header = chosen['header']
                seq = chosen['cds_seq'] or extract_sequence_at_header(cds_path, raw_header)
                if not seq:
                    nomatch_records.append(
                        f"no matches\t{cds_path.name}\t{prot_header}\tDeterministically chosen header had no sequence: {raw_header}"
                    )
                    continue
                out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
                ofh.write(out_header + "\n")
                ofh.write(wrap(seq) + "\n")
                unamb_matched += 1
                resolved_ambig += 1
                continue

            # single clear hit
            raw_header = hits[0]
            seq = extract_sequence_at_header(cds_path, raw_header)
            if not seq:
                nomatch_records.append(f"no matches\t{cds_path.name}\t{prot_header}\tMatched header but no sequence: {raw_header}")
                continue
            out_header = f">{species_key}_{raw_header[1:]}" if raw_header.startswith(">") else f">{species_key}_{raw_header}"
            ofh.write(out_header + "\n")
            ofh.write(wrap(seq) + "\n")
            unamb_matched += 1

    # write logs
    if nomatch_records:
        nomatch_path.parent.mkdir(parents=True, exist_ok=True)
        with open(nomatch_path, "w", encoding="utf-8") as nm:
            nm.write("TAG\tFILE\tPROT_HEADER\tEXTRA\n")
            for r in nomatch_records:
                nm.write(r + "\n")
    else:
        nomatch_path = None

    if ambcalls_records:
        ambcalls_path.parent.mkdir(parents=True, exist_ok=True)
        with open(ambcalls_path, "w", encoding="utf-8") as ac:
            ac.write(
                "species_key\tfile\tprotein_header\tchosen_header\tchosen_len\tothers_count\tother_headers\tother_lens\t"
                "reason\ttie_path\tquery_aa\torig_aa\tall_headers\tall_aa_lens\n"
            )
            for r in ambcalls_records:
                ac.write(r + "\n")
    else:
        ambcalls_path = None

    return unamb_matched, resolved_ambig, total, out_path, ambcalls_path, nomatch_path

def resolve_ambig_by_longest(fa_path: Path, raw_headers: List[str]) -> Tuple[str, List[str], List[int]]:
    if not raw_headers:
        raise ValueError("No headers to resolve.")
    if len(raw_headers) == 1:
        return raw_headers[0], [], []
    pairs = [(seq_len_for_header(fa_path, h), h) for h in raw_headers]
    pairs.sort(key=lambda x: (-x[0],))
    chosen = pairs[0][1]
    others = [h for (_L, h) in pairs[1:]]
    other_lens = [L for (L, _h) in pairs[1:]]
    return chosen, others, other_lens

# =========================
# CLI
# =========================
def main():
    ap = argparse.ArgumentParser(
        description="Extract CDS sequences per protein using EXACT ID matches; resolves Ensembl multi-hits via original protein files."
    )
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--proteins_fasta", help="Single protein FASTA/alignment to process.")
    g.add_argument("--in_dir", help="Directory of protein FASTAs to batch process.")
    ap.add_argument("--glob", default="*.fa*", help="Glob when using --in_dir")
    ap.add_argument("--cds_dir", required=True, help="Directory containing CDS FASTAs (.fa, .fna; may be gzipped).")
    ap.add_argument("--out_dir", help="Directory for final *_CDS.fasta outputs.")
    ap.add_argument("--ambcalls_dir", help="Directory for per-file .ambcalls.tsv logs.")
    ap.add_argument("--no_match_log_dir", help="Directory for *.nomatch.txt logs.")
    ap.add_argument("--primary-tie-break", dest="primary_tie_break", action="store_true",
                    help="If multiple exact matches, pick the longest and log to ambcalls TSV.")
    ap.add_argument("--species_json", required=True,
                    help="species_mapping.json (species_key -> { matched_cds_file, scheme{type,class}, ... })")
    ap.add_argument("--species_csv", help="(Optional) CSV with 'CDS' column; overrides CDS filenames from JSON if provided.")
    ap.add_argument("--orig_proteins_dir", help="(Optional) Directory with original Ensembl protein FASTAs; enables Ensembl transcript-based resolution.")
    ap.add_argument("--debug", action="store_true", help="Verbose debug prints of candidates tried")
    args = ap.parse_args()

    cds_dir = Path(args.cds_dir)

    # Ensure these are always defined (prevents UnboundLocalError)
    out_dir = Path(args.out_dir).resolve() if args.out_dir else None
    ambcalls_dir = Path(args.ambcalls_dir).resolve() if args.ambcalls_dir else None
    no_match_log_dir = Path(args.no_match_log_dir).resolve() if args.no_match_log_dir else None
    orig_dir = None  # <--- IMPORTANT: define even if --orig_proteins_dir is omitted

    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Writing outputs to: {out_dir}")

    if ambcalls_dir:
        ambcalls_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Ambcalls logs in: {ambcalls_dir}")

    if no_match_log_dir:
        no_match_log_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] No-match logs in: {no_match_log_dir}")

    # Load species JSON (source of truth)
    with open(args.species_json, "r", encoding="utf-8") as jh:
        sj: Dict[str, Any] = json.load(jh)

    # Build scheme lookup and species->CDS filename map
    scheme_lookup: Dict[str, Dict[str, Any]] = {}
    species_to_cds: Dict[str, str] = {}

    for skey, rec in sj.items():
        scheme_lookup[skey] = rec.get("scheme") or {"type": "ens_only"}
        cds_name = rec.get("matched_cds_file") or rec.get("expected_cds_file") or ""
        if cds_name:
            species_to_cds[skey] = Path(cds_name).name

    # Optional: override CDS filenames via CSV mapping (if provided)
    if args.species_csv:
        print("[INFO] Overriding CDS filenames with species CSV mapping")
        species_map_csv = read_species_to_cds_map(Path(args.species_csv))
        species_to_cds.update(species_map_csv)

    # Optional: original Ensembl protein directory -> build map
    if args.orig_proteins_dir:
        orig_dir = Path(args.orig_proteins_dir).resolve()
        if not orig_dir.exists():
            raise SystemExit(f"--orig_proteins_dir not found: {orig_dir}")
        build_species_to_origprot(orig_dir)
        print(f"[INFO] Ensembl transcript recovery enabled from: {orig_dir}  (indexed {len(SPECIES_TO_ORIGPROT)} species keys)")

    if args.proteins_fasta:
        unamb, ambres, total, out_path, ambcalls_path, nm_path = process_one_file(
            Path(args.proteins_fasta), cds_dir, out_dir,
            ambcalls_dir, no_match_log_dir,
            primary_tie_break=args.primary_tie_break,
            scheme_lookup=scheme_lookup,
            species_to_cds=species_to_cds,
            debug=args.debug,
            orig_proteins_dir=orig_dir,
        )
        msg = f"[{Path(args.proteins_fasta).name}] matched {unamb} (resolved {ambres} multiple hits) / {total}  ->  {out_path}"
        if ambcalls_path:
            msg += f"\n  Ambcalls TSV: {ambcalls_path}"
        if nm_path:
            msg += f"\n  No-match log: {nm_path}"
        print(msg)
        return

    in_dir = Path(args.in_dir)
    files = sorted(in_dir.glob(args.glob))
    if not files:
        raise SystemExit(f"No files matched in {in_dir} with pattern {args.glob}")
    for f in files:
        unamb, ambres, total, out_path, ambcalls_path, nm_path = process_one_file(
            f, cds_dir, out_dir,
            ambcalls_dir, no_match_log_dir,
            primary_tie_break=args.primary_tie_break,
            scheme_lookup=scheme_lookup,
            species_to_cds=species_to_cds,
            debug=args.debug,
            orig_proteins_dir=orig_dir,
        )
        msg = f"[{f.name}] matched {unamb} (resolved {ambres} multiple hits) / {total}  ->  {out_path}"
        if ambcalls_path:
            msg += f"  (ambcalls: {ambcalls_path.name})"
        if nm_path:
            msg += f"  (no-match: {nm_path.name})"
        print(msg)

if __name__ == "__main__":
    main()
