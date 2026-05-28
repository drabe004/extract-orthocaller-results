#!/bin/bash
#SBATCH --job-name=copy
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=2g
#SBATCH --cpus-per-task=1

cd ${BASE_DIR}/Orthocaller_Codon_Alignments

# ---- PATHS ----
CDS_IN="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5_CDS/FRAME_CORRECTED"

CDS_OUT="${CDS_IN}/SIFTED_CDS"

PROT_IN="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_V8_ShortestDist_NoBranchReassignments5_CDS/Translated_Proteins/FRAME_FIXED"
PROT_OUT="${PROT_IN}/SIFTED_TRANSLATED_CDS"

# ---- MAKE OUTPUT DIRS ----
mkdir -p "$CDS_OUT" "$PROT_OUT"

echo "Copying missing CDS files..."
for f in "$CDS_IN"/*.fa "$CDS_IN"/*.fna "$CDS_IN"/*.fasta "$CDS_IN"/*.fas; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    if [[ ! -e "$CDS_OUT/$base" ]]; then
        cp -p "$f" "$CDS_OUT/"
    fi
done

echo "Copying missing protein files..."
for f in "$PROT_IN"/*.fa "$PROT_IN"/*.faa "$PROT_IN"/*.fasta "$PROT_IN"/*.fas; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    if [[ ! -e "$PROT_OUT/$base" ]]; then
        cp -p "$f" "$PROT_OUT/"
    fi
done

echo "Done copying missing files."
