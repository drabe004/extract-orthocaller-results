#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --tmp=5g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=drabe004@umn.edu 




cd ${BASE_DIR}/Orthocaller_Codon_Alignments


module load conda 
source activate orthocaller



python count_species_pal2nal.py \
  -i ${BASE_DIR}/Orthocaller_Codon_Alignments/PAL2NAL_alns_V8_ShortestDist_NoBranchReassignments_SIFTED_May7_cleanedforSelTestsNoStopCodons/ \
  -c ${BASE_DIR}/Species_List_FG.txt \
  -o pal2nal_species_counts_V8_ShortestDist_NoBranchReassignments_May7.csv