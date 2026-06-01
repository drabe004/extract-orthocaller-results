#!/bin/bash -l
#SBATCH -J get_orig_prots
#SBATCH --array=1-4000
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6g
#SBATCH --tmp=4g
#SBATCH --mail-type=FAIL

###############################################################################
# Original Protein Sequence Recovery
#
# Retrieves full-length source protein sequences corresponding to
# Orthocaller-derived orthogroup FASTA files. Each SLURM array task processes
# one orthogroup, resolves species-specific protein identifiers using a genome
# key, and extracts matching primary transcript sequences from the reference
# proteome dataset.
#
# No-match outputs are written separately to support troubleshooting,
# provenance tracking, and iterative recovery of unresolved sequences.
#
# Author: Danielle Drabeck
###############################################################################

set -euo pipefail

# -------------------------------
# Logging
# -------------------------------
LOGDIR="${BASE_DIR}/logs"
mkdir -p "$LOGDIR"

exec > "${LOGDIR}/get_orig_prots_v1.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOGDIR}/get_orig_prots_v1.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# -------------------------------
# Environment
# -------------------------------
module load compatibility/agate-centos7
module load conda
source activate orthocaller

# -------------------------------
# Inputs / Outputs
# -------------------------------
FILELIST="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins/filelist.txt"
PRIMARY_ROOT="Path/To/Genomic/Proteomes"
SPECIES_KEY="${BASE_DIR}/Genome_Key.csv" ##Key file that lists all species--genomes--proteomes

ODIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins/ORIGINALSEQS"
NMDIR="${ODIR}/nomatch"

SCRIPT="${BASE_DIR}/ExtractOrthocallerResults/GetOriginalProts_V1.py"

mkdir -p "$ODIR" "$NMDIR"

# -------------------------------
# Select one orthogroup for this array task
# -------------------------------
: "${SLURM_ARRAY_TASK_ID:?This script must be run as a SLURM array task}"

IDX="$SLURM_ARRAY_TASK_ID"
FILE=$(sed -n "${IDX}p" "$FILELIST" | tr -d '\r' || true)

if [[ -z "${FILE:-}" ]]; then
  echo "[WARN] No file found for array index ${IDX}"
  exit 0
fi

if [[ ! -f "$FILE" ]]; then
  echo "[ERR] Input FASTA not found: $FILE"
  exit 2
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "[ERR] Python script not found: $SCRIPT"
  exit 2
fi

echo "[INFO] Task index: $IDX"
echo "[INFO] Input FASTA: $FILE"
echo "[INFO] Primary transcript root: $PRIMARY_ROOT"
echo "[INFO] Species key: $SPECIES_KEY"
echo "[INFO] Output directory: $ODIR"
echo "[INFO] No-match directory: $NMDIR"

# -------------------------------
# Run original protein recovery
# -------------------------------
python3 "$SCRIPT" \
  --in_fasta "$FILE" \
  --primary_dir "$PRIMARY_ROOT" \
  --species_key_csv "$SPECIES_KEY" \
  --out_dir "$ODIR" \
  --no_match_dir "$NMDIR" \
  --debug

echo "[DONE] $FILE"