#!/bin/bash -l
#SBATCH -J translate
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL


###############################################################################
# CDS Translation Pipeline
#
# Translates orthogroup-specific coding sequence (CDS) datasets into protein
# sequences using the standard genetic code. Protein translations are
# generated for all recovered CDS files and organized into a dedicated output
# directory for downstream frame validation, sequence quality control, and
# protein alignment workflows.
#
# This step provides an independent translation-based verification layer that
# enables detection of frame shifts, premature stop codons, sequence retrieval
# errors, and other inconsistencies prior to codon-aware evolutionary
# analyses.
#
# Author: Danielle Drabeck
###############################################################################





set -euo pipefail

cd ${BASE_DIR}/ExtractOrthocallerResults

module load conda
source activate orthocaller


python translate_CDS_folder.py \
  --indir ${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_CDS \
  --outdir ${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_CDS/Translated_Proteins
