#!/bin/bash
#SBATCH -J realign_mafft
#SBATCH -A mcgaughs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --array=3306,6530
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drabe004@umn.edu
#SBATCH -o ${BASE_DIR}/ExtractOrthocallerResults/MAFFT_logfiles_translatedCDSMay4/realign_%A_%a.out
#SBATCH -e ${BASE_DIR}/ExtractOrthocallerResults/MAFFT_logfiles_translatedCDSMay4/realign_%A_%a.err

set -euo pipefail
set -euo pipefail

module load mafft

# filelist has filenames only
ALIGNMENT_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5_CDS/Translated_Proteins/FRAME_FIXED/SIFTED_TRANSLATED_CDS"
ALIST="${ALIGNMENT_DIR}/protein_files.list"

# Output
BASE_OUTPUT_DIR="${ALIGNMENT_DIR}/ReAligned_MAFFT"
mkdir -p "${BASE_OUTPUT_DIR}"

# Derive the file for this array index (sed is 1-based)
ALIGNMENT_FILE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ALIST}" | tr -d '\r')"

if [[ -z "${ALIGNMENT_FILE}" ]]; then
  echo "No filename for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}. Check ${ALIST} length."
  exit 1
fi

input_file="${ALIGNMENT_DIR}/${ALIGNMENT_FILE}"

if [[ ! -s "${input_file}" ]]; then
  echo "Input file not found or empty: ${input_file}"
  exit 2
fi

# Clean base name: strip only the LAST extension
fname="$(basename "${input_file}")"
base_name="${fname%.*}"

output_file="${BASE_OUTPUT_DIR}/${base_name}_realigned.fa"

echo "[$(date)] Host: ${HOSTNAME}"
echo "[$(date)] MAFFT starting on: ${input_file}"
echo "[$(date)] Writing to: ${output_file}"
echo "CPUs requested: ${SLURM_CPUS_PER_TASK:-1}"

mafft \
  --auto \
  --anysymbol \
  --thread "${SLURM_CPUS_PER_TASK:-1}" \
  "${input_file}" \
  > "${output_file}"

echo "[$(date)] Realignment complete for ${input_file}"
