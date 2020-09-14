#!/usr/bin/env bash
# This script downloads the (manually) selected BED files for Laura's ENCODE correlation analyses.
# Alastair Droop, 2020-04-09

# Set the download BED folder:
OUT_DIR="./bed"
mkdir -p ${OUT_DIR}

# Set the input file:
IN_FILE="./encode-histone-files.txt"

# Loop through the file:
while read LINE
do
	if [[ $LINE == "#"* ]]; then continue; fi # Remove comment lines
	LINE_D=($LINE)
	OUT_FILE="${OUT_DIR}/${LINE_D[0]}-${LINE_D[1]}.bed.gz"
	aws --no-sign-request s3 cp ${LINE_D[2]} ${OUT_FILE}
	gunzip ${OUT_FILE}
done < ${IN_FILE}
