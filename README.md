## Illumina metagenomic analysis manual (IMAM).

Required files can be found here:

- workflow/Snakefile.smk: This is the heart of the workflow. 
- scripts/prepare_project.py: Use this to setup your project directory
- envs/environment.yml: Use this to recreate the necessary conda environment

To view the manual in its entirety, please visit the [EMC-Viroscience website](https://lucvzon.github.io/EMC-Viroscience.github.io/workflows.html)!

## Brief workflow summary:
1. **Quality control**
  - Deduplication: [cd-hit-dup](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki)
  - Trimming: ... [fast](https://github.com/OpenGene/fastp)
  
2. De novo assembly

3. Taxonomic classification

4. Work work work . . .


## Installation and quick start