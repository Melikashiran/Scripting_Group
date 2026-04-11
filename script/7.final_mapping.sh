#!/bin/bash


source /apps/profiles/modules_asax.sh.dyn
module load hisat2
module load samtools
module load stringtie
module load gffread

BASE_DIR="/home/aubpma001/prdm9/Scripting_Group/script"
RAW_DIR="/home/aubpma001/prdm9/Scripting_Group/script/raw_data"
MAP_DIR="/home/aubpma001/prdm9/Scripting_Group/script/alignments"
COUNT_DIR="/home/aubpma001/prdm9/Scripting_Group/script/counts"

THREADS=8

mkdir -p "$MAP_DIR" "$COUNT_DIR"

echo "Starting multi-species RNA-seq pipeline..."

for species_dir in "$RAW_DIR"/*; do
    [ -d "$species_dir" ] || continue
    species=$(basename "$species_dir")

    echo "======================================"
    echo "Species: $species"
    echo "======================================"

    ########## Locate assembly directory
ASSEMBLY_DIR="$species_dir/assembly"

GENOME_GZ=$(find "$ASSEMBLY_DIR" -maxdepth 1 -type f \( -name "*.fna.gz" -o -name "*.fa.gz" \) | head -n 1)
GFF_GZ=$(find "$ASSEMBLY_DIR" -maxdepth 1 -type f -name "*.gff.gz" | head -n 1)

if [[ -z "$GENOME_GZ" || -z "$GFF_GZ" ]]; then
    echo "Missing genome or GFF for $species → skipping"
    continue
fi

########## Define standardized outputs
REF_FASTA="$ASSEMBLY_DIR/${species}.fasta"
REF_GTF="$ASSEMBLY_DIR/${species}.gtf"
INDEX_PREFIX="$ASSEMBLY_DIR/${species}_index"
SS_FILE="$ASSEMBLY_DIR/${species}.ss"
EXON_FILE="$ASSEMBLY_DIR/${species}.exon"

########## Prepare reference (ONLY ONCE)
if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
    echo "Preparing reference for $species..."

    ########## Unzip genome (only once)
    if [[ ! -f "$REF_FASTA" ]]; then
        gunzip -c "$GENOME_GZ" > "$REF_FASTA"
    fi

    ########## ALWAYS convert GFF → GTF
    if [[ ! -f "$REF_GTF" ]]; then
        echo "Converting GFF → GTF..."
        gunzip -c "$GFF_GZ" | gffread -T -o "$REF_GTF"
    fi

    ########## Validate GTF
    echo "Validating GTF..."
    gffread -E "$REF_GTF" 2> "$ASSEMBLY_DIR/gtf_validation.log"

    ########## Extract splice sites & exons
    hisat2_extract_splice_sites.py "$REF_GTF" > "$SS_FILE"
    hisat2_extract_exons.py "$REF_GTF" > "$EXON_FILE"

    ########## Build HISAT2 index
    hisat2-build -p $THREADS \
        --ss "$SS_FILE" \
        --exon "$EXON_FILE" \
        "$REF_FASTA" \
        "$INDEX_PREFIX"

    echo "Index built for $species"
else
    echo "Index already exists for $species"
fi

    ########## Process tissues
    for tissue_dir in "$species_dir"/*; do
        [ -d "$tissue_dir" ] || continue
        tissue=$(basename "$tissue_dir")

        for r1 in "$tissue_dir"/*_1_paired.fastq.gz; do
            [ -e "$r1" ] || continue

            r2="${r1/_1_paired.fastq.gz/_2_paired.fastq.gz}"
            [ -e "$r2" ] || continue

            sample=$(basename "$r1" _1_paired.fastq.gz)

            echo "Processing: $species | $tissue | $sample"

            SAMPLE_MAP_DIR="$MAP_DIR/$species/$tissue/$sample"
            SAMPLE_COUNT_DIR="$COUNT_DIR/$species/$tissue/$sample"

            mkdir -p "$SAMPLE_MAP_DIR" "$SAMPLE_COUNT_DIR"

            cd "$SAMPLE_MAP_DIR"

            ########## ALIGNMENT
            if [[ ! -f "$sample.sorted.bam" ]]; then

                hisat2 -p $THREADS --dta --phred33 \
                    -x "$INDEX_PREFIX" \
                    -1 "$r1" \
                    -2 "$r2" \
                    -S "$sample.sam"

                samtools view -@ $THREADS -bS "$sample.sam" > "$sample.bam"

                samtools sort -@ $THREADS "$sample.bam" -o "$sample.sorted.bam"

                samtools index "$sample.sorted.bam"

                samtools flagstat "$sample.sorted.bam" > "$sample.stats.txt"

                rm "$sample.sam" "$sample.bam"

            else
                echo "Skipping alignment (exists)"
            fi

            ########## STRINGTIE
            if [[ ! -f "$SAMPLE_COUNT_DIR/$sample.gtf" ]]; then

                stringtie -p $THREADS -e -B \
                    -G "$REF_GTF" \
                    -o "$SAMPLE_COUNT_DIR/$sample.gtf" \
                    -l "$sample" \
                    "$sample.sorted.bam"

            else
                echo "Skipping quantification (exists)"
            fi

        done
    done
done

echo "Multi-species pipeline with splice-aware indexing completed!"
