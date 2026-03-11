#!/bin/bash

# ================================================
# RNA-seq SRA Downloader
# Downloads SRAs into organized species/tissue directories
# ================================================

# Path to your CSV/TSV table: columns = Organism, Tissue, SRA
# Example: "Lepidochelys_olivacea,Liver,SRR12963788"
SRA_TABLE="sra_list.csv"

# Base directory to store raw FASTQ files
RAW_DIR="raw_data"

# Check if SRA Toolkit is installed
if ! command -v fasterq-dump &> /dev/null
then
    echo "fasterq-dump not found! Please install SRA Toolkit."
    exit 1
fi

mkdir -p "$RAW_DIR"

# Read the table line by line
while IFS=',' read -r species tissue sra; do
    # Skip header if present
    [[ "$species" == "Organism" ]] && continue

    # Clean names: replace spaces with underscores
    species_dir=$(echo "$species" | sed 's/ /_/g')
    tissue_dir=$(echo "$tissue" | sed 's/ /_/g')

    # Create directories
    mkdir -p "$RAW_DIR/$species_dir/$tissue_dir"

    # Download SRA
    echo "Downloading $sra for $species_dir / $tissue_dir..."
    
    # Use prefetch first (downloads SRA to cache), then convert to FASTQ
    prefetch "$sra" && \
    fasterq-dump "$sra" -O "$RAW_DIR/$species_dir/$tissue_dir" --split-files --gzip

    # Optional: remove the .sra file to save space
    rm -f "$sra.sra"

done < "$SRA_TABLE"

echo "All downloads completed!"
