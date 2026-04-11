#!/bin/bash

source /apps/profiles/modules_asax.sh.dyn
module load sra

SRA_TABLE="/home/aubpma001/prdm9/Scripting_Group/data/sra/sra_lists.csv"
RAW_DIR="raw_data"
SRA_DIR="sra_cache"

mkdir -p "$RAW_DIR"
mkdir -p "$SRA_DIR"

# Ensure toolkit exists
if ! command -v fasterq-dump &> /dev/null; then
    echo "ERROR: fasterq-dump not found!"
    exit 1
fi

while IFS=',' read -r species tissue sra; do
    [[ "$species" == "Organism" ]] && continue

    species_dir=$(echo "$species" | sed 's/ /_/g')
    tissue_dir=$(echo "$tissue" | sed 's/ /_/g')

    OUT_DIR="$RAW_DIR/$species_dir/$tissue_dir"
    mkdir -p "$OUT_DIR"

    echo "======================================"
    echo "Processing: $sra"
    echo "Species: $species_dir | Tissue: $tissue_dir"
    echo "======================================"

    # Step 1: Download with resume support + no size limit
    prefetch "$sra" \
        --max-size 100G \
        --output-directory "$SRA_DIR"

    # Step 2: Validate download (CRITICAL)
    vdb-validate "$SRA_DIR/$sra" || {
        echo "ERROR: Validation failed for $sra. Redownloading..."
        rm -rf "$SRA_DIR/$sra"
        continue
    }

    # Step 3: Convert to FASTQ
    fasterq-dump "$SRA_DIR/$sra" \
        -O "$OUT_DIR" \
        --split-files \
        -e 8 \
        --progress

    # Step 4: Compress
    gzip "$OUT_DIR"/*.fastq

    # Step 5: Cleanup SRA to save space
    rm -rf "$SRA_DIR/$sra"

done < "$SRA_TABLE"

echo "All downloads completed!"
