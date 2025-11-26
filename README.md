## Illumina metagenomic analysis manual (IMAM).

[![Deploy Quarto Book](https://github.com/EMC-Viroscience/illumina-metagenomic-analysis-manual/actions/workflows/deploy.yml/badge.svg)](https://github.com/EMC-Viroscience/illumina-metagenomic-analysis-manual/actions/workflows/deploy.yml)

To view the fully fledged manual, please visit the [EMC-Viroscience website](https://EMC-Viroscience.github.io/workflows.html)!

Interesting files can be found here:

- **workflow/snakefile_imam.smk**: This is the heart of the workflow. 
- **scripts/prepare_project.py**: Script to setup your project directory.
- **envs/environment.yml**: The conda environment used for the workflow.
- **envs/imam_workflow.def**: The definition file used for container creation.

## Installation and quick start

1. **Download pre-built image**
  ```
  wget https://github.com/EMC-Viroscience/illumina-metagenomic-analysis-manual/releases/latest/download/imam_workflow.sif
  ```

2. **Prepare project directory**
  ```
  singularity exec \
  --bind /mnt/viro0002:/mnt/viro0002 \
  --bind $HOME:$HOME \
  --bind $PWD:$PWD \
  imam_workflow.sif \
  python /prepare_project.py \
    -p {project.folder} \
    -n {name} \
    -r {reads} \
    --ref-genome {reference} \
    --diamond-db {database} \
    -t {threads}
  ```
  
3. **Run Snakemake workflow**
  ```
  singularity exec \
  --bind /mnt/viro0002:/mnt/viro0002 \
  --bind $HOME:$HOME \
  --bind $PWD:$PWD \
  imam_workflow.sif \
  snakemake --snakefile Snakefile \
  --cores {threads}
  ```
  
## Brief workflow summary:
1.  **Quality control**
    - Merging and decompressing with zcat.
    - **Deduplication**: Remove duplicate reads from the uncompressed FASTQ files with [cd-hit-dup](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki).
    - **Quality trimming**: Perform quality and sequence adapter trimming with [fastp](https://github.com/OpenGene/fastp).
    - **Host filtering**: Remove reads that map to a host genome (e.g., human) with [bwa](https://github.com/lh3/bwa).
  
2.  **De novo assembly**
    - Perform de novo assembly of the host-filtered reads to create contigs with [SPades](https://github.com/ablab/spades).

3.  **Taxonomic classification**
    - Annotate the aggregated contigs by assigning taxonomic classifications to them based on sequence similarity to known proteins in a database using [diamond blastx](https://github.com/bbuchfink/diamond).

4.  **Extracting Viral Sequences and Analyzing Mapped Reads**
    - Extract contigs for annotated, unannotated and viral contigs.
    - Map the quality-filtered and host-filtered reads back to the assembled contigs to quantify the abundance of different contigs in each sample.
    - Extract and count mapped reads for each annotation file.
