#!/bin/bash
#SBATCH -J cln
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 04:00:00
#SBATCH --array=4001-7709
#SBATCH -o logsMay7/cln_%A_%a.out
#SBATCH -e logsMay7/cln_%A_%a.err
#SBACTH -p astyanax,msismall

cd ${BASE_DIR}/SelectionTests

module load python
mkdir -p logsMay7


###This script deletes stop codons and the rest of the header information from the sequences bc hyphy and other programs won't tolerate those things. Header backup CSVs are written as a quick reference to the header information that was truncated each species

python cln.py \
  --indir ${BASE_DIR}/Orthocaller_Codon_Alignments/PAL2NAL_alns_V8_ShortestDist_NoBranchReassignments_SIFTED_May7/ \
  --outdir ${BASE_DIR}/Orthocaller_Codon_Alignments/PAL2NAL_alns_V8_ShortestDist_NoBranchReassignments_SIFTED_May7_cleanedforSelTestsNoStopCodons/ \
  --backupdir ${BASE_DIR}/Orthocaller_Codon_Alignments/PAL2NAL_alns_V8_ShortestDist_NoBranchReassignments_SIFTED_May7_cleanedforSelTestsNoStopCodons/header_backups \
  --array \
  --delim "_" \
  --suffix "_cln" \
  --line-width 60
