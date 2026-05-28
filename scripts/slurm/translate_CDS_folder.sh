#!/bin/bash -l
#SBATCH -J translate
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drabe004@umn.edu
set -euo pipefail

cd ${BASE_DIR}/ExtractOrthocallerResults

module load compatibility/agate-centos7
module load conda
source activate orthocaller


python translate_CDS_folder.py \
  --indir ${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault_CDS \
  --outdir ${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault_CDS/Translated_Proteins
