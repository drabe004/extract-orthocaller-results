#!/bin/bash -l
#SBATCH -J rerun_cds_nomatch
#SBATCH --array=1-4000
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --tmp=4g
#SBATCH --mail-type=FAIL

###############################################################################
# CDS Retrieval Recovery Workflow
#
# Reprocesses protein families that failed initial CDS retrieval by targeting
# entries recorded in the no-match log generated during previous extraction
# attempts. Each SLURM array task processes a single failed orthogroup and
# performs an expanded CDS search using species mappings, reference proteomes,
# and genomic CDS datasets.
#
# Outputs include recovered CDS sequences, ambiguous match reports, and
# updated no-match logs to facilitate iterative dataset completion and
# quality-control assessment.
#
# Intended as a remediation step for maximizing CDS recovery rates in
# large-scale orthology and comparative genomics workflows.
#
# Author: Danielle Drabeck
###############################################################################



set -euo pipefail

module load conda
source activate orthocaller

BASE="${BASE_DIR}/ExtractOrthocallerResults"

PROTDIR="$BASE/EXTRACTED_Proteins"
RERUN_LIST="$BASE/EXTRACTED_Proteins_CDS/nomatch/nomatchlist.txt"

ODIR="$BASE/EXTRACTED_Proteins_CDS/reruns"
NOMATCHDIR="$ODIR/nomatch"
ACDIR="$ODIR/ambcalls"
LOGDIR="$ODIR/logs"

CDSDIR="Path/To/Genomic/CDS"
SPECIES_JSON="path/to/species_mapping.json"
ORIGPROTSDIR="path/to/proteomefiles" ### If used primary transcript.py (OrthoFinder accessory script) then path to primary transcript result should be used

SCRIPT="$BASE/GetCDSFromProteinsV5_1_symbol_safe.py"

mkdir -p "$ODIR" "$NOMATCHDIR" "$ACDIR" "$LOGDIR"

exec > "$LOGDIR/rerun_cds.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "$LOGDIR/rerun_cds.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "[INFO] Starting task $SLURM_ARRAY_TASK_ID"

[[ -f "$SCRIPT" ]] || { echo "[ERR] Missing script: $SCRIPT"; exit 1; }
[[ -f "$RERUN_LIST" ]] || { echo "[ERR] Missing rerun list: $RERUN_LIST"; exit 1; }

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$RERUN_LIST" | tr -d '\r')

[[ -n "$FILE" ]] || { echo "[WARN] No file for task $SLURM_ARRAY_TASK_ID"; exit 0; }

# Convert nomatch log filename ? original protein FASTA name
FILE=$(basename "$FILE")
FILE="${FILE%_CDS.fasta.nomatch.txt}"

# Build full path
FILE="$PROTDIR/$FILE"

echo "[INFO] Protein FASTA: $FILE"

[[ -f "$FILE" ]] || { echo "[ERR] Protein FASTA not found: $FILE"; exit 2; }

python3 "$SCRIPT" \
  --proteins_fasta "$FILE" \
  --cds_dir "$CDSDIR" \
  --species_json "$SPECIES_JSON" \
  --orig_proteins_dir "$ORIGPROTSDIR" \
  --out_dir "$ODIR" \
  --ambcalls_dir "$ACDIR" \
  --no_match_log_dir "$NOMATCHDIR" \
  --debug

echo "[DONE] $FILE"