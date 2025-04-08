#!/usr/bin/env python

import os, sys
import argparse
import shutil
import pandas as pd
import glob
import subprocess
import yaml

parser = argparse.ArgumentParser(description="Interactive tool for setting up a Snakemake project.")
parser.add_argument("-p", "--project_dir", help="Project directory path (default: current directory)", default=".") # current directory
parser.add_argument("-n", "--study_name", required=True, help="Name of the study")
parser.add_argument("-r", "--raw_fastq_dir", required=True, help="Directory containing raw FASTQ files")
parser.add_argument("-t", "--threads", required=False, help="Maximum number of threads for the Snakefile (default: 8)", default=8)

args = parser.parse_args() # reads command from user

# 0. Project directory path
project_dir = args.project_dir
if not os.path.exists(project_dir):
    os.makedirs(project_dir)
print(f"Using project directory: {project_dir}")

# 1. Copy Snakemake file and update STUDY_NAME
src_snakemake = "/snakefile_imam.smk" # This is adjusted for the singularity image
dest_snakemake = os.path.join(project_dir, "Snakefile")
shutil.copy2(src_snakemake, dest_snakemake)

with open(dest_snakemake, 'r') as file:
    filedata = file.read()

# Replace the STUDY_NAME variable
study_name = args.study_name
filedata = filedata.replace('STUDY_NAME = ""', f'STUDY_NAME = "{study_name}"')

with open(dest_snakemake, 'w') as file:
    file.write(filedata)

print(f"Copied and modified Snakemake file to: {dest_snakemake}")

 # 2. Set up raw FASTQ files
raw_fastq_dir = args.raw_fastq_dir

if not os.path.exists(raw_fastq_dir):
    print(f"Error: Raw FASTQ directory '{raw_fastq_dir}' does not exist.")
    sys.exit(1)

create_link = input("Create a 'raw_data' directory and symbolic links to the FASTQ files? (y/n): ")
if create_link.lower() == "y":
    raw_data_dir = os.path.join(project_dir, "raw_data")
    os.makedirs(raw_data_dir, exist_ok=True)

    # Iterate through the files in directory
    files = glob.glob(os.path.join(raw_fastq_dir,"*.fastq.gz"))
    for file in files:
       print(file)
       dest_file = os.path.join(raw_data_dir,os.path.basename(file))
       # Check if the file link already exists before creating the link
       if not os.path.exists(dest_file):
            # Make source file path absolute for the link
            abs_src_file = os.path.abspath(file)
            command=f'ln -s "{abs_src_file}" "{dest_file}"' # Use absolute source
            subprocess.run(command, shell=True, check=True)
       else:
            print(f"Link: {dest_file} already exists, skipping.")
            continue

    print(f"Created 'raw_data' directory and symbolic links to FASTQ files.")

# 3. Generate sample configuration (TSV)
raw_fastq_dir = args.raw_fastq_dir
# Ensure raw_fastq_dir itself is absolute for robust path joining and globbing
raw_fastq_dir = os.path.abspath(raw_fastq_dir)

r1_files = glob.glob(os.path.join(raw_fastq_dir, "*_R1_001.fastq.gz"))
sample_data = []
for r1_file_path in r1_files: # Use a different variable name to avoid confusion
    r1_file_abs = os.path.abspath(r1_file_path)
    # Extract sample name from the original path or the absolute one, shouldn't matter
    sample_name = os.path.splitext(os.path.basename(r1_file_path))[0].split("_R1")[0]

    # Construct the expected R2 file path relative to the (now absolute) raw_fastq_dir
    r2_file_path = os.path.join(raw_fastq_dir, f"{sample_name}_R2_001.fastq.gz")

    if not os.path.exists(r2_file_path):
        print(f"Warning: Corresponding R2 file not found for {sample_name}")
        continue

    # *** FIX HERE: Make the R2 path absolute before appending ***
    r2_file_abs = os.path.abspath(r2_file_path)

    sample_data.append({"sample": sample_name, "fq1": r1_file_abs, "fq2": r2_file_abs}) # Use absolute paths for both

df = pd.DataFrame(sample_data)
tsv_file = os.path.join(project_dir, "sample.tsv")
df.to_csv(tsv_file, sep="\t", index=False)

print(f"Created sample configuration file: {tsv_file}")

# 4. Generate general settings configuration (YAML)
ref_genome = input("Enter path to reference genome (.fna): ")
# Check that the input .fna exists in the directory
if not os.path.exists(ref_genome):
    print(f"Error: The file {ref_genome} does not exist. Please specify your input reference genome")
    sys.exit(1)
diamond_db = input("Enter path to DIAMOND database (.dmnd): ")
# Check that the diamond exists in the directory
if not os.path.exists(diamond_db):
    print(f"Error: The file {diamond_db} does not exist. Please specify your input diamond database")
    sys.exit(1)

# 5. Set threads
threads = int(args.threads)

config_data = {
    "ref_genome": os.path.abspath(ref_genome),
    "diamond_db": os.path.abspath(diamond_db),
    "threads": threads
}

yaml_file = os.path.join(project_dir, "config.yaml")
with open(yaml_file, "w") as outfile:
    yaml.dump(config_data, outfile, default_flow_style=False)

print(f"Created general settings configuration file: {yaml_file}")
