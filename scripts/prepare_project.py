#!/usr/bin/env python

import os, sys
import argparse
import shutil
import re
import pandas as pd
import glob
import subprocess
import yaml
import json # version checker
import urllib.request # version checker
from pathlib import Path # version checker

# ==================================
#           UPDATE CHECKER
# ==================================

def check_for_updates(repo_owner: str, repo_name: str):
    """
    Checks for a new version of the software on GitHub.

    This function compares the local version (from /IMAM_VERSION) with the latest
    release on GitHub. If a newer version is available, it notifies the user
    and asks if they want to abort the script.

    Args:
        repo_owner (str): The owner of the GitHub repository (e.g., 'your_username').
        repo_name (str): The name of the GitHub repository (e.g., 'naam-workflow').
    """
    # --- Define some colors for the output message for better visibility ---
    class Colors:
        YELLOW = '\033[93m'
        GREEN = '\033[92m'
        RED = '\033[91m'
        ENDC = '\033[0m'

    # Define the two possible locations for the version file.
    container_path = Path("/IMAM_VERSION")

    # When running locally, the script is in '.../scripts/amplicon_project.py'.
    # The version file is in the project root, so we go up two directories.
    script_location = Path(__file__).resolve() # Absolute path to this script
    project_root = script_location.parent.parent # .../scripts/ -> .../
    local_path = project_root / "IMAM_VERSION"

    # Now, figure out which one to use.
    if container_path.exists():
        local_version_file = container_path
    elif local_path.exists():
        local_version_file = local_path
    else:
        # If neither file exists, we can't check, so we just skip.
        print("Warning: IMAM_VERSION file not found in standard locations. Skipping update check.", file=sys.stderr)
        return

    local_version_str = local_version_file.read_text().strip()

    # 2. Fetch the latest release from the GitHub API
    api_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/releases/latest"
    try:
        with urllib.request.urlopen(api_url, timeout=5) as response:
            if response.status != 200:
                # Handle cases like private repos or API rate limiting
                return
            data = json.loads(response.read().decode())
            latest_version_str = data['tag_name']

    except Exception:
        # If there's any issue (no internet, API down), just fail silently.
        # The user's workflow should not be blocked by the version check.
        return

    # 3. Compare the versions
    # A simple helper to parse versions like 'v1.2.3' into comparable tuples (1, 2, 3)
    def parse_version(v_str):
        return tuple(map(int, v_str.lstrip('v').split('.')))

    local_v = parse_version(local_version_str)
    latest_v = parse_version(latest_version_str)

    if latest_v > local_v:
        # A newer version is available. Prompt the user.
        # --- CHANGE 1: Add ENDC to the separator lines ---
        print(f"\n{Colors.YELLOW}{'='*70}{Colors.ENDC}", file=sys.stderr)
        print(f"{Colors.YELLOW}  WARNING: A new version of this workflow is available!{Colors.ENDC}", file=sys.stderr)
        print(f"{Colors.YELLOW}{'-'*70}{Colors.ENDC}", file=sys.stderr)
        print(f"  You are using version: {Colors.RED}{local_version_str}{Colors.ENDC}", file=sys.stderr)
        print(f"  Latest version is:     {Colors.GREEN}{latest_version_str}{Colors.ENDC}", file=sys.stderr)
        print(f"  Please update at: https://github.com/{repo_owner}/{repo_name}/releases", file=sys.stderr)
        # --- Also add ENDC to the final separator line ---
        print(f"{Colors.YELLOW}{'='*70}{Colors.ENDC}", file=sys.stderr)
        
        try:
            # --- CHANGE 2: Color the input prompt string itself, ending with ENDC ---
            prompt = f"  Do you want to abort the current setup to update? (y/n): {Colors.ENDC}"
            answer = input(prompt).lower().strip()

            if answer == 'y':
                print(f"\n{Colors.RED}Aborting script. Please pull the latest container version.{Colors.ENDC}", file=sys.stderr)
                sys.exit(1)
            else:
                # --- CHANGE 3: Add ENDC to the final message for a clean exit ---
                print(f"Continuing with the current version...\n{Colors.ENDC}", file=sys.stderr)
        except (EOFError, KeyboardInterrupt):
            # Handle Ctrl+D or Ctrl+C during input as an abort
            print(f"\n\n{Colors.RED}Aborting script.{Colors.ENDC}", file=sys.stderr)
            sys.exit(1)

parser = argparse.ArgumentParser(description="Interactive tool for setting up a Snakemake project.")
parser.add_argument("-p", "--project_dir", help="Project directory path (default: current directory)", default=".") # current directory
parser.add_argument("-n", "--study_name", required=True, help="Name of the study")
parser.add_argument("-r", "--raw_fastq_dir", required=True, help="Directory containing raw FASTQ files")
parser.add_argument("-t", "--threads", type=int, required=False, help="Maximum number of threads for the Snakefile (default: 8)", default=8)
parser.add_argument("--ref-genome", required=True, help="Path to the reference genome FASTA file (.fna, .fasta, .fa)")
parser.add_argument("--diamond-db", required=True, help="Path to the DIAMOND database file (.dmnd)")

args = parser.parse_args() # reads command from user

check_for_updates(repo_owner="EMC-Viroscience", repo_name="illumina-metagenomic-analysis-manual")

# 0. Project directory path
project_dir = args.project_dir
if not os.path.exists(project_dir):
    os.makedirs(project_dir)
print(f"Using project directory: {project_dir}")

# 1. Copy Snakemake file and update STUDY_NAME
# --- Robustly find the source Snakefile ---
container_snakefile_path = Path("/snakefile_imam.smk")

# When local, it's in '<project_root>/workflow/snakefile_naam.smk'
script_location = Path(__file__).resolve()
project_root = script_location.parent.parent
local_snakefile_path = project_root / "workflow" / "snakefile_imam.smk"

if container_snakefile_path.exists():
    src_snakemake = container_snakefile_path
elif local_snakefile_path.exists():
    src_snakemake = local_snakefile_path
else:
    print(f"Error: Source Snakemake file not found at '{container_snakefile_path}' or '{local_snakefile_path}'.")
    sys.exit(1)

dest_snakemake = os.path.join(project_dir, "Snakefile")

try:
    shutil.copy2(src_snakemake, dest_snakemake)

    with open(dest_snakemake, 'r') as file:
        filedata = file.read()

    # Replace the STUDY_NAME variable
    study_name = args.study_name
    # Make the replacement more robust (e.g., handle spaces around =)
    filedata = re.sub(r'STUDY_NAME\s*=\s*""', f'STUDY_NAME = "{study_name}"', filedata)
    # filedata = filedata.replace('STUDY_NAME = ""', f'STUDY_NAME = "{study_name}"') # Original less robust way

    with open(dest_snakemake, 'w') as file:
        file.write(filedata)

    print(f"Copied and modified Snakemake file to: {dest_snakemake}")

except Exception as e:
    print(f"Error processing Snakemake file: {e}")
    sys.exit(1)

 # 2. Set up raw FASTQ files
raw_fastq_dir = args.raw_fastq_dir

if not os.path.exists(raw_fastq_dir):
    print(f"Error: Raw FASTQ directory '{raw_fastq_dir}' does not exist.")
    sys.exit(1)

# 3. Generate sample configuration (TSV)
raw_fastq_dir = args.raw_fastq_dir
# Ensure raw_fastq_dir itself is absolute for robust path joining and globbing
raw_fastq_dir = os.path.abspath(raw_fastq_dir)

r1_files = glob.glob(os.path.join(raw_fastq_dir, "*_R1_001.fastq.gz"))
sample_data = []
for r1_file_path in r1_files: # Use a different variable name to avoid confusion

    # Check if "undetermined" is in the filename
    base_filename = os.path.basename(r1_file_path)
    if "undetermined" in base_filename.lower():
        print(f" Skipping undetermined file: {base_filename}")
        continue # Skip to the next R1 file

    r1_file_abs = os.path.abspath(r1_file_path)
    # Extract sample name from the original path or the absolute one, shouldn't matter
    sample_name = os.path.splitext(os.path.basename(r1_file_path))[0].split("_R1")[0]

    # Construct the expected R2 file path relative to the (now absolute) raw_fastq_dir
    r2_file_path = os.path.join(raw_fastq_dir, f"{sample_name}_R2_001.fastq.gz")

    if not os.path.exists(r2_file_path):
        print(f"Warning: Corresponding R2 file not found for {sample_name}")
        continue

    # Make the R2 path absolute before appending
    r2_file_abs = os.path.abspath(r2_file_path)

    sample_data.append({"sample": sample_name, "fq1": r1_file_abs, "fq2": r2_file_abs}) # Use absolute paths for both

df = pd.DataFrame(sample_data)
tsv_file = os.path.join(project_dir, "sample.tsv")
df.to_csv(tsv_file, sep="\t", index=False)

print(f"Created sample configuration file: {tsv_file}")

# 4. Validate paths and generate general settings configuration (YAML)
if not os.path.exists(args.ref_genome):
    print(f"Error: The reference genome file '{args.ref_genome}' does not exist.")
    sys.exit(1)

if not os.path.exists(args.diamond_db):
    print(f"Error: The DIAMOND database file '{args.diamond_db}' does not exist.")
    sys.exit(1)

config_data = {
    "ref_genome": os.path.abspath(args.ref_genome),
    "diamond_db": os.path.abspath(args.diamond_db),
    "threads": args.threads
}

yaml_file = os.path.join(project_dir, "config.yaml")
with open(yaml_file, "w") as outfile:
    yaml.dump(config_data, outfile, default_flow_style=False)
print(f"Created general settings configuration file: {yaml_file}")

# Victory lap
print("\nProject setup complete.")
print(f"Project Directory: {project_dir}")
print(f"Snakemake File: {dest_snakemake}")
print(f"Sample Sheet: {tsv_file}")
print(f"Configuration file: {yaml_file}")
print("\nYou can now navigate to the project directory and run Snakemake.")
