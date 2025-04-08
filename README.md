## Illumina metagenomic analysis manual (IMAM).

To view the fully fledged manual, please visit the [EMC-Viroscience website](https://lucvzon.github.io/EMC-Viroscience.github.io/workflows.html)!

Interesting files can be found here:

- **workflow/Snakefile.smk**: This is the heart of the workflow. 
- **scripts/prepare_project.py**: Script to setup your project directory.
- **envs/environment.yml**: The conda environment used for the workflow.
- **envs/imam_workflow.sif**: The container image for the workflow.

## Installation and quick start

1. **Download pre-built image**
  ```
  wget https://github.com/LucvZon/illumina-metagenomic-analysis-manual/releases/tag/v1.0.0/imam_workflow.sif
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
1. **Quality control**
  - Merging and decompressing with zcat.
  - **Deduplication**: Remove duplicate reads from the uncompressed FASTQ files with [cd-hit-dup](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki){target="_blank"}.
  - **Quality trimming**: Perform quality and sequence adapter trimming with [fast](https://github.com/OpenGene/fastp){target="_blank"}.
  - **Host filtering**: Remove reads that map to a host genome (e.g., human) with [bwa](https://github.com/lh3/bwa){target="_blank"}.
  
2. **De novo assembly**
  - Perform de novo assembly of the host-filtered reads to create contigs with [SPades](https://github.com/ablab/spades){target="_blank"}.
  - Rename contigs.
  - Aggregate contigs.

3. **Taxonomic classification**
  - Annotate the aggregated contigs by assigning taxonomic classifications to them based on sequence similarity to known proteins in a database using [diamond blastx](https://github.com/bbuchfink/diamond).
  - Split the combined annotation file back into individual annotation files for each sample.
  - Parse diamond output with **post_process_diamond.py**.

4. **Extracting Viral Sequences and Analyzing Mapped Reads**
  - Extract contigs that have been annotated as viral.
  - Map the quality-filtered and host-filtered reads back to the assembled contigs to quantify the abundance of different contigs in each sample with bwa.
  - Extract mapped reads for each annotation file with samtools.
  - Count number of mapped reads with samtools.
  - Extract contigs for annotated, unannotated and viral contigs with seqkit.
