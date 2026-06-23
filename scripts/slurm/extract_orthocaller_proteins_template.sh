#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com



module load conda 
source activate orthocaller


REPO_DIR="/path/to/extract-orthocaller-results"

cd "${REPO_DIR}"

module load python


MASTER_SUMMARY="/path/to/master_summary_file.txt"
ORTHOCALLER_BASE="/path/to/Orthocaller_results"
GENERAX_KEY="/path/to/GeneRaxKey.txt"
ALN_DIR="/path/to/GeneraxFormattedAlns"
OUT_DIR="/path/to/output_directory"
  

python scripts/extract_orthocaller_proteins.py \
  --master_summary "${MASTER_SUMMARY}" \
  --orthocaller_base "${ORTHOCALLER_BASE}" \
  --generax_key "${GENERAX_KEY}" \
  --aln_dir "${ALN_DIR}" \
  --out_dir "${OUT_DIR}"