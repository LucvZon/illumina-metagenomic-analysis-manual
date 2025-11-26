# de novo assembly using metaSPades, DIAMOND taxonomic classification 


STUDY_NAME = ""

import shutil
import pandas as pd
import glob
import os
import re
import sys

# Check for an environment variable that is ONLY set inside the Singularity container.
# If it exists, we are in the container. Otherwise, we are running locally.
if "SINGULARITY_NAME" in os.environ:
    # We are in the container, scripts are in the root.
    SCRIPT_PATH = "/"
    VERSION_FILE_PATH = "/IMAM_VERSION"
    print(f"--- Running inside Singularity. Expecting version file at: {VERSION_FILE_PATH} ---", flush=True)
else:
    # We are running locally, scripts are in the 'scripts' sub-directory.
    SCRIPT_PATH = "scripts/"
    VERSION_FILE_PATH = "IMAM_VERSION"
    print(f"--- Running in a local environment. Expecting version file at: {VERSION_FILE_PATH} ---", flush=True)

# --- Print container version to the log for immediate feedback ---
try:
    with open(VERSION_FILE_PATH, "r") as f:
        container_version = f.read().strip()
    print(f"--- IMAM Workflow running with container version: {container_version} ---", flush=True)
except FileNotFoundError:
    print("--- WARNING: Could not determine container version. {VERSION_FILE_PATH} not found. ---", flush=True)

# Load general configuration
configfile: "config.yaml"

# Load sample information from the TSV file
# Ensure sample.tsv contains absolute paths in fq1, fq2 columns
try:
    # Keep the original sample name (with lane etc.) as the index for now
    sample_data = pd.read_csv("sample.tsv", sep="\t").set_index("sample", drop=False)
except Exception as e:
    print("Error reading or processing sample.tsv:")
    print(e)
    raise

# --- Reinstate logic to find BASE samples ---
def extract_base_sample_name(sample_name_with_suffix):
    """
    Extracts the base sample name, removing common sequencing suffixes
    like _S[0-9]+, _L00[1-4], _R[12], _001.
    Adjust the regex if your naming convention differs significantly.
    Example: UDI_1_S1_L001_R1_001 -> UDI_1
             SampleX_R1_001 -> SampleX
    """
    # Remove _S followed by numbers
    base_name = re.sub(r'_S\d+', '', sample_name_with_suffix)
    # Remove _L00 followed by 1-4
    base_name = re.sub(r'_L00[1-4]', '', base_name)
    # Remove _R1 or _R2
    base_name = re.sub(r'_R[12]', '', base_name)
    # Remove _001 suffix
    base_name = re.sub(r'_001$', '', base_name)
    # Remove potential leftover trailing underscores
    base_name = base_name.strip('_')
    return base_name

# Apply the function to the original 'sample' column to get base names
sample_data["base_sample"] = sample_data["sample"].apply(extract_base_sample_name)

# Get unique base sample names
BASE_SAMPLES = sorted(sample_data["base_sample"].unique().tolist())
# -----------------------------------------------

# --- Create a dictionary mapping BASE sample to its list of FASTQ files ---
#     Structure: {base_sample: {'R1': [abs_path_r1_lane1, abs_path_r1_lane2,...],
#                               'R2': [abs_path_r2_lane1, abs_path_r2_lane2,...]}}
fastq_files_for_base_sample = {}
for base_sample in BASE_SAMPLES:
    # Find all rows in the original dataframe corresponding to this base sample
    matching_rows = sample_data[sample_data["base_sample"] == base_sample]
    # Get the absolute paths from the fq1 and fq2 columns for these rows
    r1_files = sorted(matching_rows["fq1"].tolist())
    r2_files = sorted(matching_rows["fq2"].tolist())
    fastq_files_for_base_sample[base_sample] = {'R1': r1_files, 'R2': r2_files}
# -----------------------------------------------------------------------

print(f"Base Samples: {BASE_SAMPLES}")
print(f"Fastq Files Per Base Sample: {fastq_files_for_base_sample}")


rule all:
    input:
        "result/used_container_version.txt",
        # Expand based on the unique BASE samples
        expand("result/{sample}/mapping/annotated_idxstats.tsv", sample=BASE_SAMPLES),
        expand("result/{sample}/mapping/unannotated_idxstats.tsv", sample=BASE_SAMPLES),
        expand("result/{sample}/mapping/viral_idxstats.tsv", sample=BASE_SAMPLES),
        expand("result/{sample}/mapping/annotated_contigs.fasta", sample=BASE_SAMPLES),
        expand("result/{sample}/mapping/unannotated_contigs.fasta", sample=BASE_SAMPLES),
        expand("result/{sample}/mapping/viral_contigs.fasta", sample=BASE_SAMPLES),
        expand("result/readstats/{statfile}.tsv", statfile=["raw", "dedup", "qc", "filter", "mapped", "viral"]),
        "result/readstats/ALL_STATS_COMBINED.tsv",
        "result/final_processed_annotation.tsv",
        "result/final_viral_annotation.tsv"

rule record_container_version:
    message:
        "Recording the IMAM container version used for this run."
    output:
        "result/used_container_version.txt"
    shell:
        """
        cat {VERSION_FILE_PATH} > {output}
        """

rule merge_and_unzip_fastq:
    message: # Add a message for clarity
        "Merging and unzipping FASTQ files for base sample {wildcards.sample}"
    input:
        # Use the base sample name (wildcards.sample) to look up the *list* of files
        R1 = lambda wildcards: fastq_files_for_base_sample[wildcards.sample]['R1'],
        R2 = lambda wildcards: fastq_files_for_base_sample[wildcards.sample]['R2']
    output:
        # Output filenames use the base sample name
        R1 = temp("result/{sample}/{sample}_R1.fastq"),
        R2 = temp("result/{sample}/{sample}_R2.fastq")
    shell:
        """
        echo "Input R1 files for {wildcards.sample}: {input.R1}" # Debugging echo
        echo "Input R2 files for {wildcards.sample}: {input.R2}" # Debugging echo
        # zcat automatically handles multiple input files and concatenates them
        zcat {input.R1} > "{output.R1}"
        zcat {input.R2} > "{output.R2}"
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
    threads: 8
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
    threads: 4
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
    threads: 4
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
    threads: 1
    output:
        temp("result/{sample}/assembly/contigs_renamed.fasta")
    shell:
        """
        seqkit replace -p "^" -r "{wildcards.sample}_" {input} > {output}
        """
#Rename contigs with "filename"_NODE_1 etc. to be able to merge them

rule aggregate_contigs:
    input:
        expand("result/{sample}/assembly/contigs_renamed.fasta", sample = BASE_SAMPLES)
    threads: 1
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
    threads: 1
    output:
        expand("result/{sample}/annotation/diamond_output.tsv", sample = BASE_SAMPLES)
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
    threads: 1
    log: "logs/{sample}_diamond_parser.log"
    output:
        annotated = "result/{sample}/annotation/annotated_contigs.tsv",
        unannotated = "result/{sample}/annotation/unannotated_contigs.tsv"
    params:
        python_script=os.path.join(SCRIPT_PATH, "post_process_diamond_v1.0.py")
    shell:
        '''
        python {params.python_script} \
        -i {input.annotation} \
        -c {input.contigs} \
        -o {output.annotated} \
        -u {output.unannotated} \
        -log {log}
        '''

rule extract_viral_annotations:
    input:
        annotated = "result/{sample}/annotation/annotated_contigs.tsv",
    threads: 1
    output:
        viral = "result/{sample}/annotation/viral_contigs.tsv"
    shell:
        '''
        grep -E "Adnaviria|Duplodnaviria|Monodnaviria|Riboviria|Ribozyviria|Varidnaviria" {input.annotated} > {output.viral} || touch {output.viral}
        '''

rule map_reads_to_contigs:
    input:
        R1 = "result/{sample}/filtered/{sample}_R1.fastq",
        R2 = "result/{sample}/filtered/{sample}_R2.fastq",
        S = "result/{sample}/filtered/{sample}_S.fastq",
        contigs = ancient("result/{sample}/assembly/contigs.fasta")
    output:
        "result/{sample}/mapping/contigs.bam"
    threads: 4
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
    threads: 1
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
    threads: 1
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
    threads: 1
    shell:
        """
        seqkit grep -f <(cut -f1 {input.annotated}) {input.contigs} > {output.annotated}
        seqkit grep -f <(cut -f1 {input.unannotated}) {input.contigs} > {output.unannotated}
        seqkit grep -f <(cut -f1 {input.viral}) {input.contigs} > {output.viral}
        """

### READ STATISTICS ###

# --- Step 1: Count reads per sample at each major step ---

rule count_reads_raw:
    input:
        R1 = "result/{sample}/{sample}_R1.fastq",
        R2 = "result/{sample}/{sample}_R2.fastq"
    output:
        temp("result/{sample}/readstats/raw.count")
    shell:
        """
        R1_COUNT=$(( $(cat {input.R1} | wc -l) / 4 ))
        R2_COUNT=$(( $(cat {input.R2} | wc -l) / 4 ))
        echo -e "{wildcards.sample}\t$R1_COUNT\t$R2_COUNT" > {output}
        """

rule count_reads_dedup:
    input:
        R1 = "result/{sample}/dedup/{sample}_R1.fastq",
        R2 = "result/{sample}/dedup/{sample}_R2.fastq"
    output:
        temp("result/{sample}/readstats/dedup.count")
    shell:
        """
        R1_COUNT=$(( $(cat {input.R1} | wc -l) / 4 ))
        R2_COUNT=$(( $(cat {input.R2} | wc -l) / 4 ))
        echo -e "{wildcards.sample}\t$R1_COUNT\t$R2_COUNT" > {output}
        """

rule count_reads_qc:
    input:
        R1 = "result/{sample}/dedup_qc/{sample}_R1.fastq",
        R2 = "result/{sample}/dedup_qc/{sample}_R2.fastq",
        S = "result/{sample}/dedup_qc/{sample}_S.fastq"
    output:
        temp("result/{sample}/readstats/qc.count")
    shell:
        """
        R1_COUNT=$(( $(cat {input.R1} | wc -l) / 4 ))
        R2_COUNT=$(( $(cat {input.R2} | wc -l) / 4 ))
        S_COUNT=$(( $(cat {input.S} | wc -l) / 4 ))
        echo -e "{wildcards.sample}\t$R1_COUNT\t$R2_COUNT\t$S_COUNT" > {output}
        """

rule count_reads_filter:
    input:
        R1 = "result/{sample}/filtered/{sample}_R1.fastq",
        R2 = "result/{sample}/filtered/{sample}_R2.fastq",
        S = "result/{sample}/filtered/{sample}_S.fastq"
    output:
        temp("result/{sample}/readstats/filter.count")
    shell:
        """
        R1_COUNT=$(( $(cat {input.R1} | wc -l) / 4 ))
        R2_COUNT=$(( $(cat {input.R2} | wc -l) / 4 ))
        S_COUNT=$(( $(cat {input.S} | wc -l) / 4 ))
        echo -e "{wildcards.sample}\t$R1_COUNT\t$R2_COUNT\t$S_COUNT" > {output}
        """

rule count_mapped_reads:
    input:
        "result/{sample}/mapping/contigs.bam"
    output:
        temp("result/{sample}/readstats/mapped.count")
    shell:
        """
        # Count primary mapped reads only
        # -F 2308: exclude supplementary, secondary, unmapped
        # -f 65: require read 1, paired
        # -f 129: require read 2, paired
        # -F 2309: exclude supplementary, secondary, unmapped, paired (i.e. get singletons)
        R1_COUNT=$(samtools view -c -F 2308 -f 65 {input})
        R2_COUNT=$(samtools view -c -F 2308 -f 129 {input})
        S_COUNT=$(samtools view -c -F 2309 {input})
        echo -e "{wildcards.sample}\t$R1_COUNT\t$R2_COUNT\t$S_COUNT" > {output}
        """

rule count_reads_viral:
    input:
        "result/{sample}/mapping/viral_contigs.bam"
    output:
        temp("result/{sample}/readstats/viral.count")
    shell:
        """
        # Count primary mapped reads only, using the same logic as for all mapped reads
        R1_COUNT=$(samtools view -c -F 2308 -f 65 {input})
        R2_COUNT=$(samtools view -c -F 2308 -f 129 {input})
        S_COUNT=$(samtools view -c -F 2309 {input})
        echo -e "{wildcards.sample}\t$R1_COUNT\t$R2_COUNT\t$S_COUNT" > {output}
        """

# --- Step 2: Aggregate counts from all samples into per-step reports ---

rule aggregate_pe_stats:
    input:
        expand("result/{sample}/readstats/{step}.count", sample=BASE_SAMPLES, step="{step}")
    output:
        "result/readstats/{step,raw|dedup}.tsv"
    shell:
        """
        echo -e "Sample\tR1\tR2\tTotal" > {output}
        cat {input} | awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $2+$3}}' >> {output}
        """

rule aggregate_pes_stats:
    input:
        expand("result/{sample}/readstats/{step}.count", sample=BASE_SAMPLES, step="{step}")
    output:
        "result/readstats/{step,qc|filter|mapped|viral}.tsv"
    shell:
        """
        echo -e "Sample\tR1\tR2\tS\tTotal" > {output}
        cat {input} | awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $4, $2+$3+$4}}' >> {output}
        """

# --- Step 3: Combine all per-step reports into a single summary file ---

rule combine_all_stats:
    input:
        raw="result/readstats/raw.tsv",
        dedup="result/readstats/dedup.tsv",
        qc="result/readstats/qc.tsv",
        filter="result/readstats/filter.tsv",
        mapped="result/readstats/mapped.tsv",
        viral="result/readstats/viral.tsv"
    output:
        "result/readstats/ALL_STATS_COMBINED.tsv"
    shell:
        """
        echo -e "Sample\tStep\tR1\tR2\tS\tTotal" > {output}
        # Process files with 2 read types (R1, R2)
        tail -n +2 {input.raw} | awk 'BEGIN {{OFS="\\t"}} {{print $1, "raw", $2, $3, "-", $4}}' >> {output}
        tail -n +2 {input.dedup} | awk 'BEGIN {{OFS="\\t"}} {{print $1, "dedup", $2, $3, "-", $4}}' >> {output}
        # Process files with 3 read types (R1, R2, S)
        tail -n +2 {input.qc} | awk 'BEGIN {{OFS="\\t"}} {{print $1, "qc", $2, $3, $4, $5}}' >> {output}
        tail -n +2 {input.filter} | awk 'BEGIN {{OFS="\\t"}} {{print $1, "filter", $2, $3, $4, $5}}' >> {output}
        tail -n +2 {input.mapped} | awk 'BEGIN {{OFS="\\t"}} {{print $1, "mapped", $2, $3, $4, $5}}' >> {output}
        tail -n +2 {input.viral} | awk 'BEGIN {{OFS="\\t"}} {{print $1, "viral", $2, $3, $4, $5}}' >> {output}
        """

### FINAL ANNOTATION REPORT ###

rule create_per_sample_annotation_with_counts:
    input:
        annotation = "result/{sample}/annotation/annotated_contigs.tsv",
        counts = "result/{sample}/mapping/annotated_idxstats.tsv"
    output:
        temp("result/{sample}/annotation/annotated_with_counts.tsv")
    shell:
        """
        # This awk script is now robust to changes in the number of columns from the annotation file.
        awk 'BEGIN {{ FS="\\t"; OFS="\\t" }}
             FNR==NR && FNR > 1 {{ read_counts[$1] = $2; next }}
             FNR!=NR {{
                # Extract contig length from the contig ID (column 1)
                split($1, parts, "_");
                contig_len = parts[4];
                
                # Get the read count, defaulting to 0 if not found
                count = ($1 in read_counts) ? read_counts[$1] : 0;
                
                # Print the new leading columns: Sample, original Contig_ID, extracted Contig_Length and Read_Count
                printf "%s\\t%s\\t%s\\t%s", "{wildcards.sample}", $1, contig_len, count;
                
                # Loop through and print all the original annotation columns, from column 2 to the end
                for (i=2; i<=NF; i++) {{
                    printf "\\t%s", $i;
                }}
                
                # Print a newline character to finish the line
                printf "\\n";
             }}' {input.counts} {input.annotation} > {output}
        """

rule aggregate_final_annotation:
    input:
        expand("result/{sample}/annotation/annotated_with_counts.tsv", sample=BASE_SAMPLES)
    output:
        "result/final_processed_annotation.tsv"
    shell:
        """
        # Create the full, correct header reflecting all columns from DIAMOND, the python script, and our additions.
        # The order MUST match the output of the awk script in the rule above.
        echo -e "Sample\\tContig_ID\\tContig_Length\\tRead_Count\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstaxids\\ttaxon_name\\tspecies\\tgenus\\tfamily\\torder\\tclass\\tphylum\\tkingdom\\tsuperkingdom" > {output}
        
        # Concatenate all the per-sample files into the final report.
        cat {input} >> {output}
        """

rule aggregate_viral_annotation:
    input:
        "result/final_processed_annotation.tsv"
    output:
        "result/final_viral_annotation.tsv"
    shell:
        """
        echo -e "Sample\\tContig_ID\\tContig_Length\\tRead_Count\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstaxids\\ttaxon_name\\tspecies\\tgenus\\tfamily\\torder\\tclass\\tphylum\\tkingdom\\tsuperkingdom" > {output}
        
        # Concatenate all the per-sample files into the final report.
        grep -E "Adnaviria|Duplodnaviria|Monodnaviria|Riboviria|Ribozyviria|Varidnaviria" {input} >> {output} || touch {output}
        """
		
### DESTROY .snakemake/ AFTER THE WORKFLOW HAS SUCCESFULLY FINISHED
onsuccess:
    """
    This code runs only after the entire workflow completes successfully.
    It launches a detached, background process to remove the .snakemake
    directory a few seconds after the main workflow process exits.
    This avoids the "file handle in use" race condition.
    """
    import subprocess
    import shlex

    print("Workflow finished successfully.")
    print("Spawning a detached cleanup process to remove .snakemake directory in a few seconds...")

    cleanup_script = "import time, shutil; time.sleep(2); shutil.rmtree('.snakemake', ignore_errors=True)"
    command = shlex.split(f'{sys.executable} -c "{cleanup_script}"')
    
    subprocess.Popen(
        command,
        start_new_session=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    
onerror:
    """
    Optional: This code runs if the workflow fails at any point.
    Useful for explicitly stating that cleanup WON'T happen.
    """
    print("Workflow failed. The .snakemake directory will NOT be removed for debugging.")
