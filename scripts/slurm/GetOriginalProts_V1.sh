#!/bin/bash -l
#SBATCH -J get_orig_prots_v1
#SBATCH --array=1-4000
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6g
#SBATCH --tmp=4g
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=drabe004@umn.edu

set -euo pipefail

### Logging
LOGDIR="${BASE_DIR}/ExtractOrthocallerResults/logs_getoriginalprots_v1_BRAV9"
mkdir -p "$LOGDIR"
exec > "${LOGDIR}/get_orig_prots_v1.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOGDIR}/get_orig_prots_v1.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

### Environment
module load compatibility/agate-centos7
module load conda
source activate orthocaller

### Inputs / Outputs
FILELIST="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault/Unaligned_Fastas/filelist.txt"
PRIMARY_ROOT="/projects/standard/mcgaughs/drabe004/Orthofinder_Datasets/MASTER_FISH_PROTS_128sp"
SPECIES_KEY="${BASE_DIR}/125_Fish_Genome_Key.csv"

ODIR="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V9_BRAdefault/ORIGINALSEQS_Unaligned"
NMDIR="$ODIR/nomatch"
mkdir -p "$ODIR" "$NMDIR"

### One file per task
: "${SLURM_ARRAY_TASK_ID:?This script must be run as a SLURM array task}"
IDX="$SLURM_ARRAY_TASK_ID"

FILE=$(sed -n "${IDX}p" "$FILELIST" || true)
if [[ -z "${FILE:-}" ]]; then
  echo "[WARN] No file for index $IDX"
  exit 0
fi
FILE="${FILE%$'\r'}"
if [[ ! -f "$FILE" ]]; then
  echo "[ERR] Not found: '$FILE'"
  exit 2
fi

echo "[INFO] Task $IDX processing: $FILE"
echo "[INFO] PRIMARY_ROOT=$PRIMARY_ROOT (expects 'primary_transcripts/' under it)"
echo "[INFO] SPECIES_KEY=$SPECIES_KEY"
echo "[INFO] Output -> $ODIR"
ls
### Run
SCRIPT="GetOriginalProts_V1.py"

python3 "$SCRIPT" \
  --in_fasta "$FILE" \
  --primary_dir "$PRIMARY_ROOT" \
  --species_key_csv "$SPECIES_KEY" \
  --out_dir "$ODIR" \
  --no_match_dir "$NMDIR" \
  --debug

echo "[DONE] $FILE"
