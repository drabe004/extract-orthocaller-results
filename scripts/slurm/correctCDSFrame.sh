#!/bin/bash
#SBATCH --job-name=sift_frames
#SBATCH --output=${BASE_DIR}/CDSframeLOG_%A_%a.out
#SBATCH --error=${BASE_DIR}/CDSframeLOG_%A_%a.err
#SBATCH --array=1-4000
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=2g
#SBATCH --cpus-per-task=1

set -euo pipefail

module load conda
source activate orthocaller

###############################################################################
# CDS Reading Frame Correction
#
# Performs automated reading-frame correction of coding sequences using
# validated translated protein sequences as references. Each SLURM array task
# processes a single orthogroup, identifies frame discrepancies between CDS
# and protein datasets, and generates corrected CDS sequences suitable for
# downstream codon-aware alignment and evolutionary analyses.
#
# Outputs include:
#   1. Frame-corrected CDS FASTA files.
#   2. Per-orthogroup correction reports documenting all modifications.
#
# This workflow serves as a quality-control step to recover coding sequences
# affected by frame shifts, translation inconsistencies, or sequence retrieval
# artifacts prior to comparative genomics and molecular evolution analyses.
#
# Author: Danielle Drabeck
###############################################################################


BASE_DIR="${BASE_DIR}"

WORKDIR="${BASE_DIR}/Orthocaller_Codon_Alignments"
LOG_DIR="${WORKDIR}/CDSframefix"
mkdir -p "$LOG_DIR"
cd "$WORKDIR"

# -----------------------------
# INPUT DIRS / LIST
# -----------------------------
PROT_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_GetCDSResults/Translated_CDS/FRAME_CORRECTED"
PROT_LIST="${PROT_DIR}/protein_files.list"

CDS_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_GetCDSResults"

# -----------------------------
# OUTPUT DIR
# -----------------------------
OUT_DIR="${CDS_DIR}/FRAME_CORRECTED"
mkdir -p "$OUT_DIR"

# -----------------------------
# GET FILENAME FOR THIS TASK
# -----------------------------
FILENAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROT_LIST" | tr -d '\r')"

if [[ -z "$FILENAME" ]]; then
  echo "ERROR: No filename found for task ID ${SLURM_ARRAY_TASK_ID} in $PROT_LIST"
  exit 1
fi

PROT_FILE="${PROT_DIR}/${FILENAME}"
BASE="${FILENAME%.faa_protein.faa}"

CDS_FILE="${CDS_DIR}/${BASE}.faa_CDS.fasta"
OUT_CDS="${OUT_DIR}/${BASE}.faa_CDS.FRAME_CORRECTED.fasta"
REPORT="${OUT_DIR}/${BASE}.frame_correct_report.csv"

# -----------------------------
# CHECK INPUTS
# -----------------------------
if [[ ! -s "$PROT_FILE" ]]; then
  echo "MISSING PROTEIN: $PROT_FILE"
  exit 0
fi

if [[ ! -s "$CDS_FILE" ]]; then
  echo "MISSING CDS: $CDS_FILE"
  exit 0
fi

# OPTIONAL restart-safe skip:
# if [[ -s "$OUT_CDS" && -s "$REPORT" ]]; then
#   echo "SKIP (already done): $BASE"
#   exit 0
# fi

# -----------------------------
# RUN
# -----------------------------
python3 correctCDSFrame.py \
  --cds_fasta "$CDS_FILE" \
  --protein_fasta "$PROT_FILE" \
  --out_cds_fasta "$OUT_CDS" \
  --report_csv "$REPORT" \
  --genetic_code_table 1

echo "Finished $BASE"
