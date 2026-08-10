import gzip
import sys
import os
import re

def open_fastq(filepath):
    """Dynamically open compressed or uncompressed FASTQ files."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'rt')

def sanitize_filename(name):
    """Replaces spaces and invalid characters with underscores."""
    name = re.sub(r'[^A-Za-z0-9]+', '_', name)
    return name.strip('_')

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 bin_reads.py <kaiju_names.out> <fq1.gz> <fq2.gz> <output_dir>")
        sys.exit(1)

    kaiju_file = sys.argv[1]
    fq1_path = sys.argv[2]
    fq2_path = sys.argv[3]
    out_dir = sys.argv[4]

    # Target taxonomic levels and their corresponding indices in the kaiju string
    # String format: superkingdom(0), phylum(1), class(2), order(3), family(4), genus(5), species(6)
    target_levels = {
        'family': 4,
        'genus': 5,
        'species': 6
    }

    # Step 1: Create subdirectories for the taxonomic levels
    for level in target_levels.keys():
        os.makedirs(os.path.join(out_dir, level), exist_ok=True)

    # Step 2: Map Read IDs to Taxonomic levels
    print(f"Mapping reads from {os.path.basename(kaiju_file)}...")
    read_to_taxa = {}
    
    with open(kaiju_file, 'r') as f:
        for line in f:
            if line.startswith('C'):  # Only process Classified reads
                parts = line.strip('\n').split('\t')
                
                if len(parts) >= 4:
                    read_id = parts[1]
                    tax_string = parts[3]
                    
                    # Split by ';' and strip whitespace
                    tax_ranks = [rank.strip() for rank in tax_string.split(';')]
                    
                    # Ensure the string actually has enough ranks
                    if len(tax_ranks) >= 7:
                        read_to_taxa[read_id] = {}
                        
                        for level, index in target_levels.items():
                            taxon_name = tax_ranks[index]
                            if taxon_name and taxon_name != "NA":
                                safe_name = sanitize_filename(taxon_name)
                                read_to_taxa[read_id][level] = safe_name

    if not read_to_taxa:
        print("No classified reads found. Skipping binning.")
        return

    # Step 3: Parse FASTQ files and bin to FASTA
    print("Binning reads into FASTA files...")
    
    # Dictionary to hold open file handles for each level and taxon
    file_handles = { 'family': {}, 'genus': {}, 'species': {} }
    
    try:
        with open_fastq(fq1_path) as f1, open_fastq(fq2_path) as f2:
            while True:
                # Read 4 lines for R1
                h1 = f1.readline().strip()
                if not h1: break  # EOF
                s1 = f1.readline().strip()
                f1.readline(); f1.readline()  # Skip + and quality
                
                # Read 4 lines for R2
                h2 = f2.readline().strip()
                s2 = f2.readline().strip()
                f2.readline(); f2.readline()

                # Extract exact Read ID
                read_id = h1.split()[0][1:]

                if read_id in read_to_taxa:
                    # Write to all applicable taxonomic levels for this read
                    for level, taxon in read_to_taxa[read_id].items():
                        
                        # Open the file if it hasn't been opened yet
                        if taxon not in file_handles[level]:
                            fasta_path = os.path.join(out_dir, level, f"{taxon}.fasta")
                            file_handles[level][taxon] = open(fasta_path, 'a')
                        
                        # Write R1 and R2
                        file_handles[level][taxon].write(f">{h1[1:]}\n{s1}\n")
                        file_handles[level][taxon].write(f">{h2[1:]}\n{s2}\n")

    finally:
        # Ensure all file handles are closed properly
        counts = {'family': 0, 'genus': 0, 'species': 0}
        for level in file_handles:
            for f in file_handles[level].values():
                f.close()
                counts[level] += 1
                
    print(f"Created bins in {out_dir}:")
    print(f"  -> {counts['family']} Family FASTA files")
    print(f"  -> {counts['genus']} Genus FASTA files")
    print(f"  -> {counts['species']} Species FASTA files")

if __name__ == "__main__":
    main()
