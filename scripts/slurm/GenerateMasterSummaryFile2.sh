#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email
#SBATCH --job-name=master_summary

###############################################################################
# Orthocaller Master Summary Generation
#
# Aggregates orthogroup classifications from Orthocaller outputs into a
# consolidated master summary table. Orthogroups are filtered according to
# user-defined foreground and background species representation thresholds,
# enabling rapid identification of gene families suitable for downstream
# comparative genomics and evolutionary analyses.
#
# The resulting summary file serves as the primary index for subsequent
# sequence extraction, alignment generation, and selection-testing workflows.
#
# Author: Danielle Drabeck
###############################################################################



cd ${BASE_DIR}/Orthocaller

# Load environment
module load python



python GenerateMasterSummaryFile2.py \
  -i ${BASE_DIR}/Orthocaller_Danielle/Orthocaller_Output_Dir/ \
  -o master_summary_file.txt \
  --min-cavefish 3 \ ### Change to your species of choice
  --min-background 30  ## Change to your species of choice