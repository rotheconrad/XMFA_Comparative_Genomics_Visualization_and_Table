
# XMFA Comparative Genomics Visualization

![Python](https://img.shields.io/badge/python-3.11-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-active-success)

### Visualizing multi-genome alignments with gene and repeat annotations

This repository provides a Python tool for generating **publication-quality visualizations of multi-genome alignments** produced by **Mauve / progressiveMauve (XMFA format)**.

The script integrates:

* genome alignment differences
* gene models
* repeat annotations
* SNP variation across genomes

into a single multi-track figure useful for **comparative genomics and candidate locus analysis**.

The output figures are **fully vectorized (PDF/SVG)** and can be edited in **Adobe Illustrator**, **Inkscape**, or other vector graphics software.

The README provided here explains the workflow, required inputs, and how to interpret the resulting figures. 

---

# Table of Contents

* [Example Output](#example-output)
* [Quick Start](#quick-start)
* [Workflow Overview](#workflow-overview)
* [Step 1 — Extract Genomic Regions](#step-1--extract-genomic-regions)
* [Step 2 — Generate XMFA Alignment](#step-2--generate-xmfa-alignment)
* [Step 3 — Retrieve Annotation Files](#step-3--retrieve-annotation-files)
* [Step 4 — Visualization Script](#step-4--visualization-script)
* [Conda Environment](#conda-environment)
* [Example Usage](#example-usage)
* [Output Files](#output-files)
* [Visualization Tracks](#visualization-tracks)
* [Interpreting the Visualization](#interpreting-the-visualization)
* [Why Use XMFA Instead of VCF](#why-use-xmfa-instead-of-vcf)
* [Performance Notes](#performance-notes)
* [Repository Structure](#repository-structure)
* [Applications](#applications)
* [Future Improvements](#future-improvements)
* [Citation](#citation)
* [License](#license)

---

# Example Output

Below is an example visualization showing variation around a candidate promoter region.

```
[IMAGE PLACEHOLDER: docs/example_alignment.png]
```

The visualization contains four tracks:

1. **Gene track** — gene models and exon structure
2. **Repeat track** — transposable element annotations
3. **SNP track** — polymorphic alignment columns
4. **Alignment matrix** — nucleotide differences relative to the reference genome

Each row represents a genome, allowing rapid identification of **haplotypes, structural variants, and promoter polymorphisms**.

---

# Quick Start

### 1️⃣ Extract a genomic region

```bash
samtools faidx genome.fa Gm15:1-4000000 > genome_chr15_1_4Mb.fa
```

### 2️⃣ Align genomes with Mauve

```bash
progressiveMauve \
  --output=Gm15_alignment.xmfa \
  genome1.fa genome2.fa genome3.fa
```

### 3️⃣ Generate the visualization

```bash
python plot_xmfa_window.py \
  --xmfa Gm15_alignment.xmfa \
  --ref-sid Wm82a6 \
  --seqid Gm15 \
  --start 2050000 \
  --end 2080000 \
  --genes-gff genes.gff3.gz \
  --repeats-gff repeats.gff3.gz \
  --order Wm82a6 Benning Clark Harosoy Kingawa Minsoy MiyakoWhite \
  --out promoter_region.pdf
```

---

# Workflow Overview

The typical workflow used in this repository is illustrated below.

```
[IMAGE PLACEHOLDER: docs/workflow_diagram.png]
```

```
Genome assemblies (FASTA)
        │
        ▼
Extract genomic region
(samtools faidx)
        │
        ▼
Region FASTA files
(one per genome)
        │
        ▼
Multiple genome alignment
(progressiveMauve)
        │
        ▼
XMFA alignment file
        │
        ▼
Retrieve annotations
(GFF3 gene + repeat tracks)
        │
        ▼
plot_xmfa_window.py
        │
        ▼
Visualization output
(PDF + SVG + summary tables)
```

---

# Step 1 — Extract Genomic Regions

Large genome assemblies can be hundreds of megabases or gigabases long.
For faster alignment and visualization it is often helpful to extract a **smaller genomic window**.

Example genomes:

```
Wm82a6
Benning
Clark
Harosoy
Kingawa
Minsoy
MiyakoWhite
```

Index the genome:

```bash
samtools faidx genome.fa
```

Extract a region:

```bash
samtools faidx genome.fa Gm15:1-4000000 > genome_chr15_1_4Mb.fa
```

---

# Step 2 — Generate XMFA Alignment

The visualization script uses **XMFA files produced by progressiveMauve**.

Install Mauve:

```
conda install -c bioconda mauve
```

Run alignment:

```bash
progressiveMauve \
  --output=Gm15_alignment.xmfa \
  Wm82a6_chr15_1_4Mb.fa \
  Benning_chr15_1_4Mb.fa \
  Clark_chr15_1_4Mb.fa \
  Harosoy_chr15_1_4Mb.fa \
  Kingawa_chr15_1_4Mb.fa \
  Minsoy_chr15_1_4Mb.fa \
  MiyakoWhite_chr15_1_4Mb.fa
```

---

# Step 3 — Retrieve Annotation Files

To provide biological context, the visualization overlays:

* gene annotations
* exon structures
* repeat annotations

Annotation sources include:

```
Phytozome
SoyBase
Ensembl Plants
```

Example annotation files:

```
Gmax_880_Wm82.a6.v1.gene_exons.gff3.gz
Gmax_880_Wm82.a6.v1.repeatmasked_assembly_v6.0.gff3.gz
```

---

# Step 4 — Visualization Script

The script

```
plot_xmfa_window.py
```

generates a multi-track visualization consisting of:

1. Gene track
2. Repeat track
3. SNP track
4. Alignment difference matrix

---

# Conda Environment

Create environment:

```bash
mamba create -n xmfa_plot python=3.11
mamba activate xmfa_plot
```

Install required packages:

```bash
mamba install numpy pandas matplotlib
```

Optional preprocessing tools:

```
samtools
mauve
```

---

# Example Usage

Example command for plotting a promoter region:

```
python plot_xmfa_window.py \
  --xmfa Gm15_alignment.xmfa \
  --ref-sid Wm82a6 \
  --seqid Gm15 \
  --start 2050000 \
  --end 2080000 \
  --genes-gff Gmax_880_Wm82.a6.v1.gene_exons.gff3.gz \
  --repeats-gff Gmax_880_Wm82.a6.v1.repeatmasked_assembly_v6.0.gff3.gz \
  --order Wm82a6 Benning Clark Harosoy Kingawa Minsoy MiyakoWhite \
  --out Glyma15G025500_promoter.pdf
```

---

# Output Files

The script produces:

### Visualization

```
region_name.pdf
region_name.svg
```

Vector graphics suitable for editing.

### Quantitative summaries

```
region_name.sample_summary.tsv
region_name.column_summary.tsv
```

These tables summarize:

* mismatch counts
* SNP counts
* gap counts
* nucleotide composition differences

---

# Visualization Tracks

## Gene Track

Displays gene structure including exons and transcription direction.

## Repeat Track

Displays transposable elements annotated by RepeatMasker.

## SNP Track

Marks alignment columns containing polymorphisms relative to the reference genome.

## Alignment Matrix

Displays nucleotide differences relative to the reference genome.

Color scheme:

| State           | Color  |
| --------------- | ------ |
| reference match | gray   |
| A mismatch      | green  |
| C mismatch      | blue   |
| G mismatch      | orange |
| T mismatch      | red    |
| N/other         | purple |
| gap             | white  |

---

# Interpreting the Visualization

This visualization integrates genome alignment with gene and repeat annotations.

Researchers can examine:

* **haplotype blocks**
* **promoter variation**
* **structural variants**
* **repeat insertions**

Patterns in the alignment matrix reveal:

* shared ancestry between lines
* recombination breakpoints
* lineage-specific mutations
* candidate regulatory variation

---

# Why Use XMFA Instead of VCF

Many genomics pipelines rely on **reference-based variant calling** stored in VCF files.

However VCF approaches can miss:

* structural variation
* transposon insertions
* large deletions
* genome rearrangements

Assembly-to-assembly alignment using **XMFA** captures:

* SNPs
* insertions and deletions
* rearrangements
* structural differences

This makes XMFA particularly useful for **crop genomes with high repeat content and structural variation**.

---

# Performance Notes

### Alignment runtime

| Region size | Genomes | Runtime    |
| ----------- | ------- | ---------- |
| 1 Mb        | 7       | ~2–5 min   |
| 4 Mb        | 7       | ~10–20 min |
| 20 Mb       | 7       | ~1–2 hr    |

### Visualization runtime

| Alignment size | Runtime |
| -------------- | ------- |
| 10k columns    | <1 sec  |
| 100k columns   | ~2 sec  |

---

# Repository Structure

```
xmfa-genome-visualization
│
├── plot_xmfa_window.py
├── README.md
│
├── docs
│   ├── example_alignment.png
│   └── workflow_diagram.png
│
├── example_data
│   ├── example_alignment.xmfa
│   ├── genes.gff3.gz
│   └── repeats.gff3.gz
│
└── example_output
    ├── example_region.pdf
    ├── example_region.svg
    ├── example_region.sample_summary.tsv
    └── example_region.column_summary.tsv
```

---

# Applications

This workflow is useful for:

* QTL candidate gene analysis
* promoter variation studies
* haplotype comparison across breeding lines
* structural variant discovery
* crop comparative genomics

---

# Future Improvements

Possible enhancements include:

* structural variant density tracks
* repeat density panels
* interactive genome browser output
* support for HAL / Cactus alignments
* haplotype clustering tools

---

# Citation

If used in publications, please cite:

```
Darling AE, Mau B, Perna NT.
progressiveMauve: multiple genome alignment with gene gain, loss and rearrangement.
PLoS One (2010)
```

---

# License

This repository is released under an GNU General Public License v3.0.

---