
# Preparing Data Resources for the eCLIP-seq Pipeline

---

## Table of Contents

- [Genome Sequence](#genome-sequence)
- [Gencode GTF with RepeatMasker](#gencode-gtf-with-repeatmasker)
- [STAR Genome Index](#star-genome-index)
- [Expressed Genes GTF](#expressed-genes-gtf)
- [CLIPper with Expressed Genes](#clipper-with-expressed-genes)

---

## Genome Sequence

Download the primary assembly genome FASTA file for **GRCh38** from Gencode:

- **File**: `GRCh38.primary_assembly.genome.fa`
- **Download**: [Gencode Human Genome](https://www.gencodegenes.org/human/)

---

## Gencode GTF with RepeatMasker

1. **Download Gencode GTF**:  
   Obtain the comprehensive annotation GTF file from Gencode:
   - **Source**: [Gencode Human Annotation](https://www.gencodegenes.org/human/)

2. **Download RepeatMasker**:  
   Download RepeatMasker annotations from UCSC Genome Browser:
   - **Source**: [UCSC RepeatMasker](https://genome.ucsc.edu/)

3. **Merge Annotations**:  
   Combine the Gencode GTF file with RepeatMasker annotations to create a unified file. Ensure that both annotations are formatted correctly for downstream processing.

---

## STAR Genome Index

Use the **STAR** aligner to generate the genome index from the downloaded FASTA and GTF files:

1. Input files:
   - `GRCh38.primary_assembly.genome.fa`
   - Gencode GTF file
   
2. Generate the STAR genome index:

---

## Expressed Genes GTF

Extract genes expressed in **DLPFC ROSMAP RNA-seq data**:

1. Input GTF file: The Gencode GTF file from the previous step.
2. Identify expressed genes using RNA-seq expression data from the **ROSMAP DLPFC** dataset.
3. Output file: `expressed_genes.gtf`, containing only the genes expressed in the DLPFC.

---

## CLIPper with Expressed Genes

Download and configure the CLIPper tool with the **ROSMAP reference** for expressed genes: `clipper/main.tar.gz`

---

