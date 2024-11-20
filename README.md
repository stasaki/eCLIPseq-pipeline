------------------------------------------------------------------------

# eCLIP-seq Processing Pipeline

This repository provides a Snakemake-based pipeline for processing eCLIP-seq data. For more details, please refer to our research paper:\
[**The YTHDF Proteins Shape the Brain Gene Signatures of Alzheimer’s Disease**](https://www.biorxiv.org/content/10.1101/2024.10.23.619425v1)

------------------------------------------------------------------------

## Table of Contents

-   [Overview](#overview)
-   [Pipeline Steps](#pipeline-steps)
    -   [Read Quality Control](#read-quality-control)
    -   [CLIPper Analysis](#clipper-analysis)
        -   [Generate Consensus Peaks](#generate-consensus-peaks)
        -   [Count Reads in Consensus Peaks](#count-reads-in-consensus-peaks)
    -   [pureCLIP Analysis](#pureclip-analysis)
        -   [Additional Read Processing](#additional-read-processing)
        -   [Generate Consensus Peaks](#generate-consensus-peaks-pureclip)
        -   [Count Reads in Consensus Peaks](#count-reads-in-consensus-peaks-pureclip)
-   [Execution](#execution)
-   [Contributing](#contributing)
-   [License](#license)

------------------------------------------------------------------------

## Overview {#overview}

This pipeline is designed to process eCLIP-seq data across multiple steps, including read quality control, peak calling with CLIPper and pureCLIP, and generating consensus peaks.

### Prerequisites
- This Snakemake pipeline is designed to run on **Google Cloud** using the **Google Life Sciences API**. If you plan to run it locally, on a cluster, or on another cloud platform, you will need to modify the configurations to suit your environment.
- Required tools:
  -   [Snakemake](https://snakemake.readthedocs.io/)
  -   [R](https://www.r-project.org/)
  -   Appropriate container images for processing.

------------------------------------------------------------------------

## Pipeline Steps {#pipeline-steps}

### Read Quality Control {#read-quality-control}

``` bash
snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s 01_process_reads/Snakefile  \
  --preemption-default 2 \
  -j 4
```

------------------------------------------------------------------------

### CLIPper Analysis {#clipper-analysis}

#### Run CLIPper for Each Sample

``` bash
snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s 02_CLIPper/Snakefile  \
  --preemption-default 2 \
  -j 4
```

#### Generate Consensus Peaks {#generate-consensus-peaks}

``` bash
Rscript ./02_CLIPper/functions/consensus_peaks.R
```

#### Count Reads in Consensus Peaks {#count-reads-in-consensus-peaks}

``` bash
snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s 02_CLIPper/Snakefile_consensus  \
  --preemption-default 2 \
  -j 4
```

------------------------------------------------------------------------

### pureCLIP Analysis {#pureclip-analysis}

#### Additional Read Processing {#additional-read-processing}

``` bash
snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s 03_pureCLIP/Snakefile_trim  \
  --preemption-default 2 \
  -j 4
```

#### Run pureCLIP for Each Sample

``` bash
snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s 03_pureCLIP/Snakefile_pureclip  \
  --preemption-default 2 \
  -j 4
```

#### Generate Consensus Peaks (pureCLIP) {#generate-consensus-peaks-pureclip}

``` bash
Rscript ./03_pureCLIP/functions/consensus_peaks.R
```

#### Count Reads in Consensus Peaks (pureCLIP) {#count-reads-in-consensus-peaks-pureclip}

``` bash
snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s 03_pureCLIP/Snakefile_pureclip_consensus  \
  --preemption-default 2 \
  -j 4
```

------------------------------------------------------------------------

## Execution {#execution}

To execute the pipeline, ensure that you have all dependencies installed and configured. Modify paths and parameters in the respective `Snakefile` as necessary. Run each step sequentially using the commands provided in the corresponding sections.

------------------------------------------------------------------------
