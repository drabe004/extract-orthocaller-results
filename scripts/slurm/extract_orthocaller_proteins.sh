#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drabe004@umn.edu




module load compatibility/agate-centos7
module load conda 
source activate orthocaller


cd /projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository/ExtractOrthocallerResults

module load python  

python extract_orthocaller_proteins.py \
  --master_summary /projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository/Orthocaller_Danielle/master_summary_file_Version9_OG__shortest_classifications_OFFICIAL_3cf_ALL.txt \
  --orthocaller_base /projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository/Orthocaller_Danielle/OFFICIAL_Version/Version9_OG__shortest_classifications/ \
  --generax_key /projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository/Generax_Scripts_and_input/GeneRaxKey2_COMBO.txt \
  --aln_dir /projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository/Generax_Scripts_and_input/GeneraxFormattedAlns2 \
  --out_dir /projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault
