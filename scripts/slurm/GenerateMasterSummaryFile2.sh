#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drabe004@umn.edu
#SBATCH --job-name=master_summary





cd ${BASE_DIR}/Orthocaller_Danielle

# Load environment
module load compatibility/agate-centos7
module load python



python GenerateMasterSummaryFile2.py \
  -i ${BASE_DIR}/Orthocaller_Danielle/OFFICIAL_Version/Version9_OG__shortest_classifications/ \
  -o master_summary_file_Version9_OG__shortest_classifications_OFFICIAL_3cf_ALL.txt \
  --min-cavefish 3 \
  --min-background 30