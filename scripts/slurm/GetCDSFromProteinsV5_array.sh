#!/bin/bash -l
#SBATCH -J getcds_array
#SBATCH --array=1-4000
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL


###############################################################################
# CDS Retrieval from Ortholog Protein Families
#
# High-throughput sequence recovery workflow for extracting coding sequences
# (CDS) corresponding to Orthocaller-derived protein families. Protein
# identifiers are matched against genomic CDS datasets using species-specific
# mappings and reference proteomes to reconstruct orthogroup-level CDS
# datasets suitable for downstream codon-aware analyses.
#
# Processing is parallelized using SLURM job arrays, with one orthogroup
# processed per task. Ambiguous matches and unresolved sequences are logged
# separately to facilitate quality control, troubleshooting, and iterative
# recovery workflows.
#
# Outputs include:
#   1. Orthogroup-specific CDS FASTA files.
#   2. Ambiguous match reports.
#   3. No-match diagnostic logs.
#
# This workflow serves as the primary CDS extraction stage linking orthology
# classifications to codon-based comparative genomics and molecular evolution
# analyses.
#
# Author: Danielle Drabeck
###############################################################################

set -euo pipefail

# -------------------------------
# Logging
# -------------------------------
LOGDIR="${BASE_DIR}/ExtractOrthocallerResults/logs"
mkdir -p "$LOGDIR"
exec > "${LOGDIR}/cds_array_v5_May12.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOGDIR}/cds_array_v5_May12.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# -------------------------------
# Environment
# -------------------------------
module load conda
source activate orthocaller

# -------------------------------
# Inputs / Outputs
# -------------------------------
FILELIST="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins/filelist.txt"

# CDS directory (CDS FASTAs used for ID lookups)
CDSDIR="path/to/genomic/CDS"

# species mapping JSON
SPECIES_JSON="Path/To/species_mapping.json"

# Original Ensembl protein FASTAs (for resolving Ensembl multi-hits by transcript)
ORIGPROTSDIR="Path/To/Genomic/Proteome/files"

# If FILELIST contains relative paths, resolve them against this dir
SHORTIN_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins"

# Output roots
ODIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_CDS"
ACDIR="$ODIR/ambcalls"
NMDIR="$ODIR/nomatch"
mkdir -p "$ODIR" "$ACDIR" "$NMDIR"

# -------------------------------
# One file per array task
# -------------------------------
: "${SLURM_ARRAY_TASK_ID:?This script must be run as a SLURM array task}"
IDX="$SLURM_ARRAY_TASK_ID"

FILE=$(sed -n "${IDX}p" "$FILELIST" || true)
if [[ -z "${FILE:-}" ]]; then
  echo "[WARN] No file for index $IDX"
  exit 0
fi

# strip potential CR (Windows line endings)
FILE="${FILE%$'\r'}"

# If entry is relative, prepend SHORTIN_DIR
if [[ ! -f "$FILE" ]]; then
  base="$(basename "$FILE")"
  alt="$SHORTIN_DIR/$base"
  if [[ -f "$alt" ]]; then
    FILE="$alt"
  else
    echo "[ERR] Not found: '$FILE' or '$alt'"
    exit 2
  fi
fi

echo "[INFO] Task $IDX processing: $FILE"
echo "[INFO] Using CDSDIR=$CDSDIR"
echo "[INFO] Using ORIGPROTSDIR=$ORIGPROTSDIR"
echo "[INFO] Writing to ODIR=$ODIR (ambcalls -> $ACDIR, nomatch -> $NMDIR)"

# -------------------------------
# Call the updated script
# (rename SCRIPT if you saved the drop-in under a different filename)
# -------------------------------
SCRIPT="GetCDSFromProteinsV5_1.py"

python3 "$SCRIPT" \
  --proteins_fasta "$FILE" \
  --cds_dir "$CDSDIR" \
  --species_json "$SPECIES_JSON" \
  --orig_proteins_dir "$ORIGPROTSDIR" \
  --out_dir "$ODIR" \
  --ambcalls_dir "$ACDIR" \
  --no_match_log_dir "$NMDIR" \
  --debug

echo "[DONE] $FILE"
