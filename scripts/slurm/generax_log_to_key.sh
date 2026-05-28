#!/bin/bash -l        
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drabe004@umn.edu
#SBATCH -p astyanax,msismall

cd /panfs/jay/groups/26/mcgaughs/drabe004/BIGFISHGENOME_DataRespository

# Set the base directory for the GeneRax run folders
BASEDIR="${BASE_DIR}/Generax_Scripts_and_input/GeneRax_Run4_Results"
OUTPUT="${BASEDIR}/GeneRaxKey4.txt"

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

echo "GeneRaxKey4.txt written to $OUTPUT"