#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --tmp=5g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=youremail.com 

###############################################################################
# Alignment Species Composition Summary
#
# Quantifies taxonomic representation across PAL2NAL-generated codon
# alignments by counting foreground and background species present in each
# orthogroup. Results are compiled into a summary table for downstream
# quality-control filtering, dataset characterization, and selection-test
# eligibility assessment.
#
# This workflow is typically used to verify species coverage requirements
# prior to comparative genomics and molecular evolution analyses.
#
# Author: Danielle Drabeck
###############################################################################


cd ${BASE_DIR}/Orthocaller_Codon_Alignments


module load conda 
source activate orthocaller



python count_species_pal2nal.py \
  -i ${BASE_DIR}/Orthocaller_Codon_Alignments/PAL2NAL_alns_SIFTED_cleanedforSelTestsNoStopCodons/ \
  -c ${BASE_DIR}/Species_List_FG.txt \
  -o pal2nal_species_counts_Version_Date.csv