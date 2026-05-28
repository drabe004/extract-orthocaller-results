#!/bin/bash -l
#SBATCH -J getcds_array_v5
#SBATCH --array=5201-7000
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drabe004@umn.edu
# (no -o/-e here; we’ll redirect at runtime)

set -euo pipefail

# -------------------------------
# Logging
# -------------------------------
LOGDIR="${BASE_DIR}/ExtractOrthocallerResults/logs_may12"
mkdir -p "$LOGDIR"
exec > "${LOGDIR}/cds_array_v5_May12.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOGDIR}/cds_array_v5_May12.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# -------------------------------
# Environment
# -------------------------------
module load compatibility/agate-centos7
module load conda
source activate orthocaller

# -------------------------------
# Inputs / Outputs
# -------------------------------
FILELIST="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault/Unaligned_Fastas/filelist.txt"

# CDS directory (CDS FASTAs used for ID lookups)
CDSDIR="/projects/standard/mcgaughs/drabe004/Orthofinder_Datasets/125_Species_OFFICIALDATASET/CDS"

# species mapping JSON
SPECIES_JSON="/projects/standard/mcgaughs/drabe004/Orthofinder_Datasets/125_Species_OFFICIALDATASET/species_mapping.json"

# Original Ensembl protein FASTAs (for resolving Ensembl multi-hits by transcript)
ORIGPROTSDIR="/projects/standard/mcgaughs/drabe004/Orthofinder_Datasets/MASTER_FISH_PROTS_128sp"

# If FILELIST contains relative paths, resolve them against this dir
SHORTIN_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault/Unaligned_Fastas"

# Output roots
ODIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault_CDS"
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
