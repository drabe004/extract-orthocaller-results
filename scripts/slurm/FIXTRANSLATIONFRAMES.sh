#!/bin/bash
#SBATCH --job-name=framefix_array
#SBATCH --output=${BASE_DIR}/ExtractOrthocallerResults/FrameFixLog2/framefix_array_%A_%a.out
#SBATCH --error=${BASE_DIR}/ExtractOrthocallerResults/FrameFixLog2/framefix_array_%A_%a.err
#SBATCH --array=4001-7950
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --tmp=4g
#SBATCH -p astyanax

set -euo pipefail

#### remind future-you where this is running
cd ${BASE_DIR}/ExtractOrthocallerResults || exit 1


#### activate environment
module load conda
source activate orthocaller

#### paths
TRANSLATED_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5_CDS/Translated_Proteins"
FRAMEFAIL_DIR="${BASE_DIR}/Orthocaller_Codon_Alignments/FrameCheckOutput_April_30/"
CDS_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5_CDS"
ORIGINAL_PROT_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5/ORIGINALSEQS_Unaligned"

OUT_DIR="${TRANSLATED_DIR}/FRAME_FIXED"
mkdir -p "${OUT_DIR}"

#### build filelist once (basenames only) -- SAFE for huge file counts
FILELIST="${TRANSLATED_DIR}/framefix_filelist.txt"
if [ ! -s "${FILELIST}" ]; then
  find "${TRANSLATED_DIR}" -maxdepth 1 -type f -name "*.faa_protein.faa" -printf "%f\n" \
    | sort > "${FILELIST}"
fi

#### SLURM array IDs are 0-based; sed is 1-based
LINE=$((SLURM_ARRAY_TASK_ID + 1))
BASENAME=$(sed -n "${LINE}p" "${FILELIST}" || true)

if [ -z "${BASENAME}" ]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} (line ${LINE}) out of range"
  exit 0
fi

INFILE="${TRANSLATED_DIR}/${BASENAME}"
echo "Processing: ${INFILE}"

#### run one file per task
python3 FIXTRANSLATIONFRAMES.py \
  --translated_dir "${TRANSLATED_DIR}" \
  --framefail_dir "${FRAMEFAIL_DIR}" \
  --cds_dir "${CDS_DIR}" \
  --original_protein_dir "${ORIGINAL_PROT_DIR}" \
  --out_dir "${OUT_DIR}" \
  --threshold 0.90 \
  --infile "${INFILE}"
