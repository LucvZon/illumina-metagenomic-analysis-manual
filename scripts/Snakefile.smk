# de novo assembly using metaSPades, DIAMOND taxonomic classification 


STUDY_NAME = ""


import pandas as pd
import glob
import os
import re

# Load general configuration
configfile: "config.yaml"

# Load sample information from the TSV file
sample_data = pd.read_csv("sample.tsv", sep="\t")

# Extract sample names without lane information
def extract_base_sample_name(sample_name):
    """Extracts the base sample name, removing lane information."""

    # Use regular expression to remove lane information (L001, L002, L003, or L004)
    base_name = re.sub(r'_L00[1-4]', '', sample_name) # Remove _L001, _L002, _L003, _L004

    # Remove _R1_001 and _R2_001
    base_name = base_name.replace("_R1_001", "").replace("_R2_001", "")
    return base_name

sample_data["base_sample"] = sample_data["sample"].apply(extract_base_sample_name)

# Get unique base sample names
base_samples = sample_data["base_sample"].unique().tolist()

# Function to find all R1 and R2 files for a given base sample
def get_all_fastq_files(base_sample):
    """Finds all R1 and R2 fastq.gz files associated with a given base sample."""
    r1_files_lane = sorted(glob.glob(f"*/{base_sample}*_L00[1-4]*_R1_001.fastq.gz"))  # With lane info
    r1_files_no_lane = sorted(glob.glob(f"*/{base_sample}*_R1_001.fastq.gz"))       # Without lane info

    r2_files_lane = sorted(glob.glob(f"*/{base_sample}*_L00[1-4]*_R2_001.fastq.gz"))  # With lane info
    r2_files_no_lane = sorted(glob.glob(f"*/{base_sample}*_R2_001.fastq.gz"))       # Without lane info


    # Combine the lists: if a file exists with a lane, use it; otherwise, use the no-lane version
    r1_files = r1_files_lane if r1_files_lane else r1_files_no_lane
    r2_files = r2_files_lane if r2_files_lane else r2_files_no_lane
    
    return r1_files, r2_files

# Create a dictionary mapping base sample to a list of R1 and R2 files
fastq_files = {}
for base_sample in base_samples:
    fastq_files[base_sample] = get_all_fastq_files(base_sample)

print(f"Base Samples: {base_samples}")
print(f"Fastq Files: {fastq_files}")


rule all:
    input:
        expand("result/{sample}/mapping/annotated_idxstats.tsv", sample = base_samples),
        expand("result/{sample}/mapping/unannotated_idxstats.tsv", sample = base_samples),
        expand("result/{sample}/mapping/viral_idxstats.tsv", sample = base_samples),
        expand("result/{sample}/mapping/annotated_contigs.fasta", sample = base_samples),
        expand("result/{sample}/mapping/unannotated_contigs.fasta", sample = base_samples),
        expand("result/{sample}/mapping/viral_contigs.fasta", sample = base_samples)

rule merge_and_unzip_fastq:
    input: 
        R1 = lambda wildcards: fastq_files[wildcards.sample][0],  # List of R1 files
        R2 = lambda wildcards: fastq_files[wildcards.sample][1]   # List of R2 files
    output:
        R1 = temp("result/{sample}/{sample}_R1.fastq"),
        R2 = temp("result/{sample}/{sample}_R2.fastq")
    shell:
        """
        zcat {input.R1} > {output.R1}
        zcat {input.R2} > {output.R2}
        """

#count the number of reads in sample expr $(cat file.fastq | wc -l) / 4
#CD-HIT-dup identifies duplicates from single or paired Illumina reads
#best to first deduplicate before QC (else you can have trimmed sequences of different length)
rule dedup_raw:
    input:
        R1 = "result/{sample}/{sample}_R1.fastq",
        R2 = "result/{sample}/{sample}_R2.fastq"
    output:
        R1 = "result/{sample}/dedup/{sample}_R1.fastq",
        R2 = "result/{sample}/dedup/{sample}_R2.fastq"
    log: "logs/{sample}_cd-hit.log"
    threads: config["threads"]
    shell:
        """
        cd-hit-dup -u 50 -i {input.R1} -i2 {input.R2} -o {output.R1} -o2 {output.R2} > {log} 2>&1
        rm result/{wildcards.sample}/dedup/*.clstr
        """

rule QC_after_dedup:
    input:
        R1 = "result/{sample}/dedup/{sample}_R1.fastq",
        R2 = "result/{sample}/dedup/{sample}_R2.fastq"
    output:
        R1 = "result/{sample}/dedup_qc/{sample}_R1.fastq",
        R2 = "result/{sample}/dedup_qc/{sample}_R2.fastq",
        S = "result/{sample}/dedup_qc/{sample}_S.fastq",
        failed = "result/{sample}/dedup_qc/{sample}_fail.fastq"
    log: "logs/{sample}_fastp.log"
    threads: config["threads"]
    shell:
        """
        fastp -i {input.R1} -I {input.R2} \
        -o {output.R1} -O {output.R2} \
        --unpaired1 {output.S} --unpaired2 {output.S} --failed_out {output.failed} \
        --length_required 50 \
        --low_complexity_filter \
        --cut_right \
        --cut_right_window_size 5 \
        --cut_right_mean_quality 25 \
        --thread {threads} \
        -j result/{wildcards.sample}/dedup_qc/qc_report.json -h result/{wildcards.sample}/dedup_qc/qc_report.html > {log} 2>&1
        """
        
#for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it.
#so file S contains reads. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])

rule filter_host:
    input:
        R1 = "result/{sample}/dedup_qc/{sample}_R1.fastq",
        R2 = "result/{sample}/dedup_qc/{sample}_R2.fastq",
        S = "result/{sample}/dedup_qc/{sample}_S.fastq"
    output:
        R1 = "result/{sample}/filtered/{sample}_R1.fastq",
        R2 = "result/{sample}/filtered/{sample}_R2.fastq",
        S = "result/{sample}/filtered/{sample}_S.fastq"
    threads: config["threads"]
    params:
        reference = config["ref_genome"]
    shell:
        """
        bwa mem -aY -t {threads} {params.reference} {input.R1} {input.R2} | \
        samtools fastq -f 4 -s /dev/null -1 {output.R1} -2 {output.R2} -
        bwa mem -aY -t {threads} {params.reference}  {input.S} | \
        samtools fastq -f 4 - > {output.S}
        """

# -aY . Y is use softclipping -a is output all alignments
# the output of bwa mem is a sam file that is piped directly into samtools in order to not save the memory consuming sam files. The samtools fastq option -c means the level of compression when writing .gz fastq files. The output bam/sam file contains information about mapped and unmapped reads. 
# -f 4 : then you obtain the unmapped (which are not host reads) reads and keep these. 
#-s option: If a singleton file is specified using the -s option then only paired sequences will be output for categories 1 and 2; paired meaning that for a given QNAME there are sequences for both category 1 and 2. If there is a sequence for only one of categories 1 or 2 then it will be diverted into the specified singletons file. This can be used to prepare fastq files for programs that cannot handle a mixture of paired and singleton reads.

rule assemble_filtered:
    input:
        R1 = "result/{sample}/filtered/{sample}_R1.fastq",
        R2 = "result/{sample}/filtered/{sample}_R2.fastq",
        S = "result/{sample}/filtered/{sample}_S.fastq"
    output:
        "result/{sample}/assembly/contigs.fasta"
    params:
        spades_folder = "result/{sample}/assembly"
    log: "logs/{sample}_spades.log"
    threads: config["threads"]
    shell:
        """
        COUNT=$(cat {input.S} | wc -l)
        if [[ $COUNT -gt 0 ]]
        then
            spades.py -t {threads} \
            --meta \
            -o {params.spades_folder} \
            -1 {input.R1} \
            -2 {input.R2} \
            -s {input.S} > {log} 2>&1
        else
            spades.py -t {threads} \
            --meta \
            -o {params.spades_folder} \
            -1 {input.R1} \
            -2 {input.R2} > {log} 2>&1
        fi        
        """
#de novo assembly of contigs 
#file with forward paired-end reads + file with reverse paired end reads
# -s file with unpaired reads 

rule rename_contigs:
    input:
        "result/{sample}/assembly/contigs.fasta"
    threads: config["threads"]
    output:
        temp("result/{sample}/assembly/contigs_renamed.fasta")
    shell:
        """
        seqkit replace -p "^" -r "{wildcards.sample}_" {input} > {output}
        """
#Rename contigs with "filename"_NODE_1 etc. to be able to merge them

rule aggregate_contigs:
    input:
        expand("result/{sample}/assembly/contigs_renamed.fasta", sample = base_samples)
    threads: config["threads"]
    output:
        "result/all_contigs.fasta"
    shell:
        """
        cat {input} > {output}
        """
#Aggregate contigs before blasting

rule annotate_contigs:
    input:
        "result/all_contigs.fasta"
    threads: config["threads"]
    params:
        database = config["diamond_db"]
    log: "logs/diamond.log"
    output:
        "result/all_contigs_annotated.tsv"
    shell:
        """
        diamond blastx \
        --frameshift 15 \
        -q {input} \
        -d {params.database} \
        -o {output} \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
        --threads {threads} \
        -b 10 -c 1 > {log} 2>&1
        """
#Annotate contigs with DIAMOND blastx

rule split_annotation_files:
    input:
        "result/all_contigs_annotated.tsv"
    threads: config["threads"]
    output:
        expand("result/{sample}/annotation/diamond_output.tsv", sample = base_samples)
    shell:
        '''
        mkdir -p tmp_split

        sed 's/_NODE/|NODE/' {input} | awk -F'|' '{{
            identifier = $1;  # Construct the identifier using the first two fields

            output_file = "tmp_split/" identifier;  # Construct the output filename

            if (!seen[identifier]++) {{
                close(output_file);  # Close the previous file (if any)
                output = output_file;  # Update the current output file
            }}

            print $2 > output;  # Append the line to the appropriate output file
        }}'

        for file in tmp_split/*; do
            mkdir -p result/$(basename ${{file}})/annotation/;
            mv $file result/$(basename ${{file}})/annotation/diamond_output.tsv;
        done

        rmdir tmp_split
        '''

rule parse_diamond_output:
    input:
        annotation = "result/{sample}/annotation/diamond_output.tsv",
        contigs = "result/{sample}/assembly/contigs.fasta"
    threads: config["threads"]
    log: "logs/{sample}_diamond_parser.log"
    output:
        annotated = "result/{sample}/annotation/annotated_contigs.tsv",
        unannotated = "result/{sample}/annotation/unannotated_contigs.tsv"
    shell:
        '''
        python /mnt/viro0002/workgroups_projects/Bioinformatics/scripts/post_process_diamond.py \
        -i {input.annotation} \
        -c {input.contigs} \
        -o {output.annotated} \
        -u {output.unannotated} \
        -log {log}
        '''

rule extract_viral_annotations:
    input:
        annotated = "result/{sample}/annotation/annotated_contigs.tsv",
    threads: config["threads"]
    output:
        viral = "result/{sample}/annotation/viral_contigs.tsv"
    shell:
        '''
        grep "Viruses$" {input.annotated} > {output.viral} || touch {output.viral}
        '''

rule map_reads_to_contigs:
    input:
        R1 = "result/{sample}/filtered/{sample}_R1.fastq",
        R2 = "result/{sample}/filtered/{sample}_R2.fastq",
        S = "result/{sample}/filtered/{sample}_S.fastq",
        contigs = ancient("result/{sample}/assembly/contigs.fasta")
    output:
        "result/{sample}/mapping/contigs.bam"
    threads: config["threads"]
    shell:
        """
        bwa index {input.contigs}
        bwa mem -Y -t {threads} {input.contigs} {input.R1} {input.R2} | samtools sort - > result/{wildcards.sample}/mapping/tmp_paired.bam
        bwa mem -Y -t {threads} {input.contigs} {input.S} | samtools sort - > result/{wildcards.sample}/mapping/tmp_singlets.bam
        samtools merge {output} result/{wildcards.sample}/mapping/tmp_paired.bam result/{wildcards.sample}/mapping/tmp_singlets.bam
        rm result/{wildcards.sample}/mapping/tmp_paired.bam result/{wildcards.sample}/mapping/tmp_singlets.bam
        """

rule extract_specific_mapped:
    input:
        annotated = "result/{sample}/annotation/annotated_contigs.tsv",
        unannotated = "result/{sample}/annotation/unannotated_contigs.tsv",
        viral = "result/{sample}/annotation/viral_contigs.tsv",
        mapped = "result/{sample}/mapping/contigs.bam"
    output:
        annotated = "result/{sample}/mapping/annotated_contigs.bam",
        unannotated = "result/{sample}/mapping/unannotated_contigs.bam",
        viral = "result/{sample}/mapping/viral_contigs.bam"
    threads: config["threads"]
    shell:
        """
        #Create temporary BED file to extract the annotated mappings from the BAM of all mappings
        cut -f1 {input.annotated} | awk -F'_' '{{print $0 "\t" 0 "\t" $4}}' > result/{wildcards.sample}/mapping/tmp.bed
        samtools view -bL result/{wildcards.sample}/mapping/tmp.bed {input.mapped} > {output.annotated}
        
        #Do the same for the unannotated contigs
        cut -f1 {input.unannotated} | awk -F'_' '{{print $0 "\t" 0 "\t" $4}}' > result/{wildcards.sample}/mapping/tmp.bed
        samtools view -bL result/{wildcards.sample}/mapping/tmp.bed {input.mapped} > {output.unannotated}

        #Do the same for the viral contigs
        cut -f1 {input.viral} | awk -F'_' '{{print $0 "\t" 0 "\t" $4}}' > result/{wildcards.sample}/mapping/tmp.bed
        samtools view -bL result/{wildcards.sample}/mapping/tmp.bed {input.mapped} > {output.viral}

        rm result/{wildcards.sample}/mapping/tmp.bed
        """

rule count_reads_mapped:
    input:
        annotated = "result/{sample}/mapping/annotated_contigs.bam",
        unannotated = "result/{sample}/mapping/unannotated_contigs.bam",
        viral = "result/{sample}/mapping/viral_contigs.bam"
    output:
        annotated = "result/{sample}/mapping/annotated_idxstats.tsv",
        unannotated = "result/{sample}/mapping/unannotated_idxstats.tsv",
        viral = "result/{sample}/mapping/viral_idxstats.tsv"
    threads: config["threads"]
    shell:
        """
        #Extract readstats, remove lines with 0 mapped reads
        samtools view -bF2052 {input.annotated} | seqkit bam -Qc - | awk '$2 != 0 {{print}}' > {output.annotated}
        samtools view -bF2052 {input.unannotated} | seqkit bam -Qc - | awk '$2 != 0 {{print}}' > {output.unannotated}
        samtools view -bF2052 {input.viral} | seqkit bam -Qc - | awk '$2 != 0 {{print}}' > {output.viral}
        """
        
rule extract_specific_contigs:
    input:
        annotated = "result/{sample}/annotation/annotated_contigs.tsv",
        unannotated = "result/{sample}/annotation/unannotated_contigs.tsv",
        viral = "result/{sample}/annotation/viral_contigs.tsv",
        contigs = "result/{sample}/assembly/contigs.fasta"
    output:
        annotated = "result/{sample}/mapping/annotated_contigs.fasta",
        unannotated = "result/{sample}/mapping/unannotated_contigs.fasta",
        viral = "result/{sample}/mapping/viral_contigs.fasta"
    threads: config["threads"]
    shell:
        """
        seqkit grep -f <(cut -f1 {input.annotated}) {input.contigs} > {output.annotated}
        seqkit grep -f <(cut -f1 {input.unannotated}) {input.contigs} > {output.unannotated}
        seqkit grep -f <(cut -f1 {input.viral}) {input.contigs} > {output.viral}
        """
