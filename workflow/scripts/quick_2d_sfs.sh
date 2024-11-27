#!/bin/bash

# Directories
vcf_dir="/nas/longleaf/home/adaigle/work/ghist_2024_work/recontact_fixed/recontact/vcfs"
output_dir="/nas/longleaf/home/adaigle/work/ghist_2024_work/recontact_fixed/recontact/sfs"
script_path="/nas/longleaf/home/adaigle/ghist_2024/workflow/scripts/compute_2d_sfs_island_mainland.py"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Export variables to use them within GNU Parallel
export output_dir
export script_path

# Function to process each file
process_file() {
    vcf_file=$1
    # Extract the replicate number from the filename (assumes "replicate_{rep}.vcf.gz" format)
    filename=$(basename "$vcf_file")
    rep=$(echo "$filename" | grep -oP '(?<=replicate_)\d+')

    # Define output file with the same replicate number
    output_file="$output_dir/replicate_${rep}_2d_sfs.tsv"

    # Run the Python script
    echo "Processing replicate $rep..."
    python "$script_path" --vcf "$vcf_file" --output "$output_file"
}

export -f process_file

# Find all VCF files and process them in parallel using 24 cores
find "$vcf_dir" -name "*.vcf.gz" | parallel -j 24 process_file
