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

###############################################################################
# Alignment Header and Stop Codon Cleanup
#
# Preprocesses codon alignments for downstream evolutionary analyses by:
#   1. Removing terminal stop codons from coding sequences.
#   2. Simplifying FASTA headers to retain only species identifiers.
#
# Many phylogenetic and molecular evolution tools (e.g., HyPhy, PAML, IQ-TREE)
# require standardized sequence names and may fail when presented with complex
# headers or embedded stop codons. Original header information is preserved in
# per-alignment backup CSV files for traceability and downstream reference.
#
# Designed for high-throughput processing using SLURM job arrays.
#
# Author: Danielle Drabeck
###############################################################################

python cln.py \
  --indir ${BASE_DIR}/Sifted_Codon_Alignments/ \
  --outdir ${BASE_DIR}/Cleaned_Codon_Alignments/ \
  --backupdir ${BASE_DIR}/Cleaned_Codon_Alns/header_backups \
  --array \
  --delim "_" \
  --suffix "_cln" \
  --line-width 60
