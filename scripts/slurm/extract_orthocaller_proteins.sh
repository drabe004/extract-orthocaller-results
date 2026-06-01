#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail

###############################################################################
# Orthocaller Sequence Extraction
#
# Extracts orthogroup-specific protein sequence datasets from Orthocaller
# classification results by integrating orthology assignments, GeneRax gene
# family mappings, and source alignments. Sequences associated with accepted
# ortholog groups are recovered and organized into analysis-ready FASTA files
# for downstream alignment, quality-control, and evolutionary analyses.
#
# This workflow serves as the primary interface between Orthocaller
# classifications and downstream comparative genomics pipelines.
#
# Author: Danielle Drabeck
###############################################################################

module load conda 
source activate orthocaller


cd ${BASE_DIR}/ExtractOrthocallerResults

module load python  

python extract_orthocaller_proteins.py \
  --master_summary ${BASE_DIR}/Orthocaller/master_summary_file.txt \
  --orthocaller_base ${BASE_DIR}/Orthocaller_classifications/ \
  --generax_key ${BASE_DIR}/Generax_Scripts_and_input/GeneRaxKey.txt \
  --aln_dir ${BASE_DIR}/Generax_Scripts_and_input/GeneraxFormattedAlns \
  --out_dir ${BASE_DIR}/ExtractOrthocallerResults/Extracted_Proteins_Version_Date
