# RNA-seq Pipeline

A modular, reproducible pipeline for multi-species RNA-seq analysis — from raw reads to differential expression.

---

##  Repository Structure
The directory structure is shown below indicting the purpose of each script
```
├── 1.download_script.sh                # Step 1 — Download raw FASTQ files
├── 2.download_assembly.sh            # Step 2 — Download reference genome & annotation
├── 3.quality_control.sh              # Step 3 — Quality control on raw reads
├── 4.trimmo_qualitycontrol.sh        # Step 4 — Adapter trimming + post-trim QC
├── 5.final_mapping.sh                # Step 6 — Final/refined mapping
├── 6.prepDE.py3                      # Step 7 — Prepare count matrix
├── 7.multi_species_rnaseq_compare.py # Step 8 — Cross-species comparison
├── metadata.csv                   # Sample metadata table
├── sra_list.csv                   # SRA table for download from NCBI repository
├── quality_control/               # Output: raw QC reports
└── trim_quality_control/          # Output: post-trim QC reports
```

---

## Dependencies
The following tools were used in the pipeline
| Tool | Purpose |
|------|---------|
| `FastQC 0.12.1` | Read quality control |
| `multiqc` | Bundles all the Fastqc html into one report |
| `Trimmomatic 0.39` | Adapter and quality trimming |
| `HISAT2 2.2.0` | Read alignment |
| `SAMtools` | BAM file handling |
| `StringTie` | Transcript assembly & quantification |
| `Python 3.x` | Downstream analysis |
| `pandas 1.5.0, matplotlib 3.6.0, numpy>=1.23.0` | Statistical analysis |

---

##  Pipeline Steps

### Step 1 — `1.download_script.sh`
**Purpose:** Downloads raw FASTQ sequencing files from a remote source (e.g., SRA, FTP, or cloud storage).

**Output:**
- Raw `.fastq.gz` files in the input data directory

---

### Step 2 — `2.download_assembly.sh`
**Purpose:** Downloads the reference genome (FASTA) and gene annotation file (GTF/GFF) needed for alignment.


**Output:**
- Reference genome `.fa` / `.fasta`
- Annotation file `.gtf` / `.gff`

---

### Step 3 — `3.quality_control.sh`
**Purpose:** Runs FastQC on raw reads to assess base quality, adapter content, and sequencing artifacts. Results are stored in `quality_control/`.


**Output:**
- FastQC HTML reports → `quality_control/`

---

### Step 4 — `4.trimmo_qualitycontrol.sh`
**Purpose:** Trims low-quality bases and adapter sequences using Trimmomatic, then runs a second round of FastQC on trimmed reads.


**Output:**
- Trimmed `.fastq.gz` files
- Post-trim QC reports → `trim_quality_control/`

---

### Step 5 — `5.final_mapping.sh`
**Purpose:** Aligns trimmed reads to the reference genome. Produces sorted BAM files.


**Output:**
- Sorted and indexed `.bam` files -> `./alignment`
/
---

### Step 6 — `6.prepDE.py3`
**Purpose:** Prepares the gene count matrix from StringTie output using the `prepDE.py` script. Generates a count table for downstream differential expression analysis.

**Input:**
- StringTie output GTF files listed in a sample list

**Output:**
- `gene_count_matrix.csv`
- `transcript_count_matrix.csv`

---

### Step  7— `7.multi_species_rnaseq_compare.py`
**Purpose:** Compares RNA-seq expression data across multiple species. Handles gene ID mapping and cross-species normalization.

**Input:**
- Count matrices per species
- `metadata.csv`

**Output:**
- Merged/normalized expression table
- DE results tables (`.csv`)
- Plots (`.png` / `.pdf`)

---

## Running the Full Pipeline

```bash
bash 1.data_download.sh
bash 2.download_assembly.sh
bash 3.quality_control.sh
bash 4.trimmo_qualitycontrol.sh
bash 6.mapping.sh
bash 7.final_mapping.sh
python3 8.prepDE.py3
python3 9.multi_species_rnaseq_compare.py
python3 10.analysis.py
```

---
[CLICK HERE TO GO BACK](https://github.com/Melikashiran/Scripting_Group/blob/main/README.md)

---
## Authors

> Melika Ghasemi Siran, Prince Mensah Ansah, Sean Onileowo, Surma Mohiudden Meem 
