#!/bin/bash

find scripts -type f \( -name "*.sh" -o -name "*.py" \) | while read file; do

    sed -i \
        's|/projects/standard/mcgaughs/drabe004/BIGFISHGENOME_DataRespository|${BASE_DIR}|g' \
        "$file"

done
