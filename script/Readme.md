# RNA-seq Pipeline

A modular, reproducible pipeline for multi-species RNA-seq analysis — from raw reads to differential expression.

---

## 📁 Repository Structure

```
├── 1.data_download.sh                # Step 1 — Download raw FASTQ files
├── 2.download_assembly.sh            # Step 2 — Download reference genome & annotation
├── 3.quality_control.sh              # Step 3 — Quality control on raw reads
├── 4.trimmo_qualitycontrol.sh        # Step 4 — Adapter trimming + post-trim QC
├── 5.download_script.sh              # Step ? — (confirm purpose)
├── 6.mapping.sh                      # Step 5 — Map reads to reference genome
├── 7.final_mapping.sh                # Step 6 — Final/refined mapping
├── 8.prepDE.py3                      # Step 7 — Prepare count matrix
├── 9.multi_species_rnaseq_compare.py # Step 8 — Cross-species comparison
├── 10.analysis.py                    # Step 9 — Downstream statistical analysis
├── metadata.csv                   # Sample metadata table
├── quality_control/               # Output: raw QC reports
└── trim_quality_control/          # Output: post-trim QC reports
```

---

## ⚙️ Dependencies

> ⚠️ *Confirm exact tools and versions from your scripts*

| Tool | Purpose |
|------|---------|
| FastQC | Read quality control |
| Trimmomatic | Adapter and quality trimming |
| HISAT2 / STAR | Read alignment |
| SAMtools | BAM file handling |
| StringTie | Transcript assembly & quantification |
| Python 3.x | Downstream analysis |
| pandas, DESeq2, etc. | Statistical analysis |

---

## 🚀 Pipeline Steps

### Step 1 — `data_download.sh`
**Purpose:** Downloads raw FASTQ sequencing files from a remote source (e.g., SRA, FTP, or cloud storage).

**Usage:**
```bash
bash data_download.sh
```

**Output:**
- Raw `.fastq.gz` files in the input data directory

---

### Step 2 — `download_assembly.sh`
**Purpose:** Downloads the reference genome (FASTA) and gene annotation file (GTF/GFF) needed for alignment.

**Usage:**
```bash
bash download_assembly.sh
```

**Output:**
- Reference genome `.fa` / `.fasta`
- Annotation file `.gtf` / `.gff`

---

### Step 3 — `quality_control.sh`
**Purpose:** Runs FastQC on raw reads to assess base quality, adapter content, and sequencing artifacts. Results are stored in `quality_control/`.

**Usage:**
```bash
bash quality_control.sh
```

**Output:**
- FastQC HTML reports → `quality_control/`

---

### Step 4 — `trimmo_qualitycontrol.sh`
**Purpose:** Trims low-quality bases and adapter sequences using Trimmomatic, then runs a second round of FastQC on trimmed reads.

**Usage:**
```bash
bash trimmo_qualitycontrol.sh
```

**Output:**
- Trimmed `.fastq.gz` files
- Post-trim QC reports → `trim_quality_control/`

---

### Step 5 — `mapping.sh`
**Purpose:** Aligns trimmed reads to the reference genome. Produces sorted BAM files.

**Usage:**
```bash
bash mapping.sh
```

**Output:**
- Sorted and indexed `.bam` files → `sample_run/`

---

### Step 6 — `final_mapping.sh`
**Purpose:** Performs the final/refined mapping step. May include filtering, deduplication, or re-alignment with updated parameters.

**Usage:**
```bash
bash final_mapping.sh
```

**Output:**
- Final `.bam` files ready for quantification

---

### Step 7 — `03b_prepDE.py3`
**Purpose:** Prepares the gene count matrix from StringTie output using the `prepDE.py` script. Generates a count table for downstream differential expression analysis.

**Usage:**
```bash
python3 03b_prepDE.py3
```

**Input:**
- StringTie output GTF files listed in a sample list

**Output:**
- `gene_count_matrix.csv`
- `transcript_count_matrix.csv`

---

### Step 8 — `multi_species_rnaseq_compare.py`
**Purpose:** Compares RNA-seq expression data across multiple species. Handles gene ID mapping and cross-species normalization.

**Usage:**
```bash
python3 multi_species_rnaseq_compare.py
```

**Input:**
- Count matrices per species
- `metadata.csv`

**Output:**
- Merged/normalized expression table

---

### Step 9 — `analysis.py`
**Purpose:** Performs downstream statistical analysis including differential expression, visualization (PCA, heatmaps, volcano plots), and result export.

**Usage:**
```bash
python3 analysis.py
```

**Input:**
- Count matrix
- `metadata.csv`

**Output:**
- DE results tables (`.csv`)
- Plots (`.png` / `.pdf`)

---

## 📋 `metadata.csv`
**Purpose:** Describes each sample — species, condition, replicate number, file paths, etc. Required by several downstream scripts.

**Expected columns** *(confirm from your file)*:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `species` | Species name |
| `condition` | Experimental condition (e.g., treatment vs control) |
| `replicate` | Replicate number |
| `file_path` | Path to FASTQ file |

---

## 🔁 Running the Full Pipeline

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

## 📝 Notes

- All paths should be configured at the top of each script before running.
- Ensure all dependencies are installed and available in your `$PATH`.
- `metadata.csv` must be updated to reflect your actual sample names and conditions.

---

## 👤 Author

> Add your name, institution, and contact here.

---

## 📄 License

> Add your license here (e.g., MIT, GPL-3.0).
