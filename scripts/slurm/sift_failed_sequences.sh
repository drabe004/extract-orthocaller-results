#!/bin/bash
#SBATCH --job-name=sift_frames
#SBATCH --output=${BASE_DIR}/logs/sift_%A_%a.out
#SBATCH --error=${BASE_DIR}/logs/sift_%A_%a.err
#SBATCH --array=1-4000
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=2g
#SBATCH --cpus-per-task=1


## PURPOSE:
## Remove sequences that still fail frame checks from CDS and translated CDS FASTAs.

## INPUTS:
## - *.frameFAIL.csv files (one per gene/OG)
## - Matching CDS FASTAs
## - Matching protein (translated CDS) FASTAs

## LOGIC:
## - For each CSV, remove `original_id` from CDS and `translated_id` from proteins
## - Match FASTA files by filename stem
## - Write cleaned FASTAs to new output directories

## MODES:
## - Array mode (--csv_index): process ONE CSV ? write ONE CDS + ONE protein FASTA
## - Full mode (no --csv_index): process all CSVs + copy unchanged FASTAs

## SAFETY:
## - Originals are never modified
## - Outputs are written only if missing (no overwrite)




set -euo pipefail

cd ${BASE_DIR}/Orthocaller_Codon_Alignments || exit 1

mkdir -p SIFT

module load conda
source activate orthocaller

FAIL_CSV_DIR="${BASE_DIR}/Orthocaller_Codon_Alignments/FrameCheckOutput"

CDS_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_CDS/FRAME_CORRECTED"

PROTEIN_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_CDS/Translated_Proteins/FRAME_CORRECTED"

OUT_CDS_DIR="${CDS_DIR}/SIFTED_CDS"
OUT_PROTEIN_DIR="${PROTEIN_DIR}/SIFTED_TRANSLATED_CDS"

python3 sift_failed_sequences.py \
  --fail_csv_dir "${FAIL_CSV_DIR}" \
  --cds_dir "${CDS_DIR}" \
  --protein_dir "${PROTEIN_DIR}" \
  --out_cds_dir "${OUT_CDS_DIR}" \
  --out_protein_dir "${OUT_PROTEIN_DIR}" \
  --pattern "*.frameFAIL.csv" \
  --remove_if_passes_is "NO" \
  --sort_csvs \
  --csv_index "${SLURM_ARRAY_TASK_ID}"
