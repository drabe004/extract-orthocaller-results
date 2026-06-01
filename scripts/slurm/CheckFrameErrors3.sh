#!/bin/bash -l
#SBATCH --job-name=frame_check
#SBATCH --output=logs/frame_check_%A_%a.out
#SBATCH --error=logs/frame_check_%A_%a.err
#SBATCH --array=1-4000
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=2g
#SBATCH --cpus-per-task=1

set -euo pipefail

###############################################################################
# Frame Check SLURM Array
#
# High-throughput quality control workflow for validating translated CDS
# sequences against reference protein alignments prior to codon-based
# evolutionary analyses. Executes one orthogroup per SLURM array task and
# identifies potential frame disruptions, translation errors, and sequence
# mismatches using a configurable similarity threshold.
#
# Intended for use in large-scale comparative genomics and phylogenomic
# pipelines.
#
# Author: Danielle Drabeck
###############################################################################

############################
# User configuration
############################

BASE_DIR="/path/to/MyBaseDir"

EXTRACTED_DIR="${BASE_DIR}/MyExtractDir"
CODON_DIR="${BASE_DIR}/MyCodonAlnDir"

PROTEIN_DATASET="MyOrthocallerExtractedProteinsDir"
CDS_DATASET="${PROTEIN_DATASET}_CDS"

FRAME_FIXED_DIR="FRAME_FIXED"
FRAMECHECK_OUTPUT="FrameCheckOutput_DateTime"

TRANSLATED_DIR="${EXTRACTED_DIR}/${CDS_DATASET}/Translated_Proteins/${FRAME_FIXED_DIR}"
ORIGINAL_DIR="${EXTRACTED_DIR}/${PROTEIN_DATASET}/ORIGINALSEQS_Unaligned"
OUT_DIR="${CODON_DIR}/${FRAMECHECK_OUTPUT}"

############################
# Setup
############################

cd "${CODON_DIR}" || exit 1

mkdir -p logs
mkdir -p "${OUT_DIR}"

module load conda
source activate orthocaller

############################
# Select file for this array task
############################

FILE=$(find "${TRANSLATED_DIR}" -maxdepth 1 -type f -name "*.faa_protein.faa" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "${FILE}" ]]; then
    echo "No translated protein file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Running frame check on:"
echo "${FILE}"

############################
# Run frame check
############################

python3 "${CODON_DIR}/CheckFrameErrors3.py" \
  --translated_dir "${TRANSLATED_DIR}" \
  --original_dir "${ORIGINAL_DIR}" \
  --out_dir "${OUT_DIR}" \
  --threshold 0.90 \
  --infile "${FILE}"