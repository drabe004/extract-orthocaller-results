#!/bin/bash
#SBATCH --job-name=copy
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=2g
#SBATCH --cpus-per-task=1

cd ${BASE_DIR}/Orthocaller_Codon_Alignments


###############################################################################
# Recover Missing Sifted Sequence Files
#
# Synchronizes CDS and translated protein datasets by copying sequence files
# that are present in the source directories but absent from the corresponding
# sifted output directories. Existing files are preserved and are not
# overwritten.
#
# This utility is primarily intended for recovery, reruns, and workflow
# maintenance following interrupted jobs or partial file generation during
# large-scale comparative genomics analyses.
#
# File metadata (timestamps and permissions) are retained during copying
# using the cp -p option.
#
# Author: Danielle Drabeck
###############################################################################



# ---- PATHS ----
CDS_IN="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_GetCDSResults/FRAME_CORRECTED"

CDS_OUT="${CDS_IN}/SIFTED_CDS"

PROT_IN="${BASE_DIR}/ExtractOrthocallerResults/EXTRACTED_Proteins_GetCDSResults/Translated_CDS/FRAME_CORRECTED"
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
