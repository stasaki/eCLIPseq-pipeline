snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s Snakefile  \
  --preemption-default 2 \
  -j 4
  
  
Rscript ./functions/consensus_peaks.R

snakemake --google-lifesciences \
  --default-remote-prefix gcbucket \
  --google-lifesciences-region us-central1 \
  --container-image snakemake/snakemake:v7.32.4 \
  --use-conda \
  --conda-frontend mamba \
  -s Snakefile_consensus  \
  --preemption-default 2 \
  -j 4
  
  