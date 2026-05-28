#!/bin/bash
#SBATCH -J cln
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 04:00:00
#SBATCH --array=1-4000
#SBATCH -o logs/cln_%A_%a.out
#SBATCH -e logs/cln_%A_%a.err

cd ${BASE_DIR}/

module load python
mkdir -p logs

###This script deletes stop codons and the rest of the header information from the sequences bc hyphy and other programs won't tolerate those things. Header backup CSVs are written as a quick reference to the header information that was truncated each species

python cln.py \
  --indir ${BASE_DIR}/Sifted_Codon_Alignments/ \
  --outdir ${BASE_DIR}/Cleaned_Codon_Alignments/ \
  --backupdir ${BASE_DIR}/Cleaned_Codon_Alns/header_backups \
  --array \
  --delim "_" \
  --suffix "_cln" \
  --line-width 60
