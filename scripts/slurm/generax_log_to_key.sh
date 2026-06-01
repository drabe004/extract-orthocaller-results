#!/bin/bash -l        
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail



###############################################################################
# GeneRax Alignment-to-Run Key Generation
#
# Builds a lookup table linking input alignment filenames to their
# corresponding GeneRax run directories by parsing GeneRax log files across
# completed analysis folders. This key file enables downstream workflows to
# connect original gene-family alignments with GeneRax output directories,
# reconciled trees, and Orthocaller classification results.
#
# Intended as a metadata indexing step for large-scale phylogenomic workflows
# where thousands of independent GeneRax runs must be tracked reproducibly.
#
# Author: Danielle Drabeck
###############################################################################


cd Your/Main/Dir

# Set the base directory for the GeneRax run folders
BASEDIR="${BASE_DIR}/GeneRax_Results"
OUTPUT="${BASEDIR}/GeneRaxKey.txt"

# Empty the output file if it exists
echo -n "" > "$OUTPUT"

# Loop through each generax.log in subdirectories
find "$BASEDIR" -name "generax.log" | while read logfile; do
    # Extract the orthogroup directory name (e.g., 1_generax)
    og_dir=$(basename "$(dirname "$logfile")")

    # Extract the alignment file from the log
    aln=$(grep -m1 "generax --families" "$logfile" | sed -E 's/.*--families[[:space:]]+([^[:space:]]+).*/\1/' | xargs basename | sed 's/\.family$/.fasta/')

    if [ -n "$aln" ]; then
        echo "$aln,$og_dir" >> "$OUTPUT"
    fi

done

echo "GeneRaxKey.txt written to $OUTPUT"