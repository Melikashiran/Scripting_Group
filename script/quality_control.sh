#!/bin/bash

#set -euo pipefail

source /apps/profiles/modules_asax.sh.dyn
module load fastqc/0.12.1
module load multiqc

RAW_DIR="/home/aubpma001/prdm9/Scripting_Group/script/raw_data"
QC_DIR="/home/aubpma001/prdm9/Scripting_Group/script/quality_control"

mkdir -p "$QC_DIR"

echo "Scanning tissue directories for FASTQ files..."

for species_dir in "$RAW_DIR"/*; do
    [ -d "$species_dir" ] || continue

    for tissue_dir in "$species_dir"/*; do
        [ -d "$tissue_dir" ] || continue

        echo "Processing: $tissue_dir"

        ########## Find FASTQ files safely
        fastq_files=$(find "$tissue_dir" -maxdepth 1 -name "*.fastq.gz")

        # Skip if none found
        [ -z "$fastq_files" ] && continue

        ########## Run FastQC in parallel
        echo "$fastq_files" | xargs -n 1 -P 8 fastqc --outdir="$tissue_dir"

    done
done

echo "Running MultiQC..."

multiqc "$RAW_DIR" -o "$QC_DIR"

tar -czvf "$QC_DIR/multiqc_report.tar.gz" -C "$QC_DIR" .

echo "QC completed successfully!"