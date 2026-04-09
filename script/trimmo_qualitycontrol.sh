#!/bin/bash

source /apps/profiles/modules_asax.sh.dyn
module load fastqc/0.12.1
module load multiqc
module load trimmomatic

RAW_DIR="/home/aubpma001/prdm9/Scripting_Group/script/raw_data"
QC_DIR="/home/aubpma001/prdm9/Scripting_Group/script/trim_quality_control"

ADAPTERS="$TRIMMOMATIC/adapters/TruSeq3-PE.fa"

mkdir -p "$QC_DIR"

echo "Starting trimming + post-trim QC..."

for species_dir in "$RAW_DIR"/*; do
    [ -d "$species_dir" ] || continue

    for tissue_dir in "$species_dir"/*; do
        [ -d "$tissue_dir" ] || continue

        echo "Processing: $tissue_dir"

        for r1 in "$tissue_dir"/*_1.fastq.gz; do
            [ -e "$r1" ] || continue

            r2="${r1/_1.fastq.gz/_2.fastq.gz}"
            [ -e "$r2" ] || continue

            sample=$(basename "$r1" _1.fastq.gz)

            echo "Trimming sample: $sample"

            out_p1="$tissue_dir/${sample}_1_paired.fastq.gz"
            out_u1="$tissue_dir/${sample}_1_unpaired.fastq.gz"
            out_p2="$tissue_dir/${sample}_2_paired.fastq.gz"
            out_u2="$tissue_dir/${sample}_2_unpaired.fastq.gz"

            ########## Skip if already trimmed
            if [[ -f "$out_p1" && -f "$out_p2" ]]; then
                echo "Skipping $sample (already trimmed)"
            else
                ########## TRIMMOMATIC
                trimmomatic PE -threads 8 \
                    "$r1" "$r2" \
                    "$out_p1" "$out_u1" \
                    "$out_p2" "$out_u2" \
                    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
                    HEADCROP:15 \
                    LEADING:3 \
                    TRAILING:3 \
                    SLIDINGWINDOW:4:15 \
                    MINLEN:36
            fi

            ########## FASTQC (ONLY on trimmed paired reads)
            fastqc "$out_p1" "$out_p2" \
                --outdir="$tissue_dir" \
                --threads 8

        done
    done
done

echo "Running MultiQC..."

multiqc "$RAW_DIR" -o "$QC_DIR"

tar -czvf "$QC_DIR/multiqc_report.tar.gz" -C "$QC_DIR" .

echo "Trimming + post-QC completed!"