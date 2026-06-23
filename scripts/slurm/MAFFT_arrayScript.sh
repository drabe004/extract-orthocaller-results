#!/bin/bash
#SBATCH -J realign_mafft
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --array=1-4000
#SBATCH --mail-type=ALL
#SBATCH -o logs/realign_%A_%a.out
#SBATCH -e logs/realign_%A_%a.err

###############################################################################
# Protein Alignment Reconstruction with MAFFT
#
# Performs de novo multiple sequence alignment of translated protein datasets
# using MAFFT. Each SLURM array task processes a single orthogroup and
# generates a standardized protein alignment for downstream codon alignment
# construction, quality-control assessment, and molecular evolution analyses.
#
# Alignments are generated using MAFFT's automatic algorithm selection mode
# (--auto), allowing the software to choose an appropriate strategy based on
# dataset size and complexity.
#
# This workflow is typically used following sequence recovery, frame
# correction, or dataset filtering steps to ensure that downstream codon
# alignments are built from consistently aligned protein sequences.
#
# Author: Danielle Drabeck
###############################################################################





set -euo pipefail
set -euo pipefail

module load mafft

# filelist has filenames only
ALIGNMENT_DIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_CDS/Translated_Proteins/FRAME_FIXED/SIFTED_TRANSLATED_CDS"
ALIST="${ALIGNMENT_DIR}/protein_files.list"

# Output
BASE_OUTPUT_DIR="${ALIGNMENT_DIR}/Realigned_MAFFT"
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
