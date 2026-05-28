#!/bin/bash -l
#SBATCH --job-name=frame_check
#SBATCH --output=${BASE_DIR}/Orthocaller_Codon_Alignments/FClogs4/framecheck_%A_%a.out
#SBATCH --error=${BASE_DIR}/Orthocaller_Codon_Alignments/FClogs4/framecheck_%A_%a.err
#SBATCH --array=4001-7950
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=2g
#SBATCH --cpus-per-task=1

set -euo pipefail

cd ${BASE_DIR}/Orthocaller_Codon_Alignments || exit 1

mkdir -p FClogs4

module load conda
source activate orthocaller

TRANSLATED_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5_CDS/Translated_Proteins/FRAME_FIXED"

ORIGINAL_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5/ORIGINALSEQS_Unaligned/"

OUT_DIR="${BASE_DIR}/Orthocaller_Codon_Alignments/FrameCheckOutput_May3"

mkdir -p "${OUT_DIR}"

FILE=$(find "${TRANSLATED_DIR}" -maxdepth 1 -type f -name "*.faa_protein.faa" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "${FILE}" ]]; then
    echo "No translated protein file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Running frame check on:"
echo "${FILE}"

python3 CheckFrameErrors3.py \
  --translated_dir "${TRANSLATED_DIR}" \
  --original_dir "${ORIGINAL_DIR}" \
  --out_dir "${OUT_DIR}" \
  --threshold 0.90 \
  --infile "${FILE}"