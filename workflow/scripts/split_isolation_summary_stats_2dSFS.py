import argparse
import numpy as np
import pandas as pd
import allel
import sys

def compute_2d_sfs(vcf_path, output_path):
    print(f"Loading VCF file from {vcf_path}...")

    try:
        # Load genotype data, variant positions, and sample names from VCF
        callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'samples'])
        genotypes = allel.GenotypeArray(callset['calldata/GT'])  # Shape: (variants, samples, ploidy)
        sample_names = callset['samples']
    
        # Separate indices for east and west samples
        east_indices = [i for i, name in enumerate(sample_names) if name.startswith('east')]
        west_indices = [i for i, name in enumerate(sample_names) if name.startswith('west')]
    
        # Genotypes for each population
        genotypes_east = genotypes[:, east_indices]
        genotypes_west = genotypes[:, west_indices]
    
    except Exception as e:
        print(f"Error reading VCF file: {e}")
        sys.exit(1)

    print("Computing 2D Site Frequency Spectrum (SFS) for east and west populations...")

    # Compute derived allele counts for each population
    ac_east = genotypes_east.count_alleles()
    ac_west = genotypes_west.count_alleles()
    dac_east = ac_east[:, 1]  # Derived allele counts for east
    dac_west = ac_west[:, 1]  # Derived allele counts for west

    # Compute 2D SFS
    sfs_2d = allel.joint_sfs(dac_east, dac_west)

    # Flatten 2D SFS and prepare for output
    sfs_flat = {}
    for i in range(sfs_2d.shape[0]):
        for j in range(sfs_2d.shape[1]):
            sfs_flat[f'sfs_{i}_{j}'] = sfs_2d[i, j]

    # Convert to DataFrame and save to output
    df = pd.DataFrame([sfs_flat])
    print(f"Saving 2D SFS to {output_path}...")
    try:
        df.to_csv(output_path, index=False, sep='\t')
    except Exception as e:
        print(f"Error saving 2D SFS: {e}")
        sys.exit(1)

    print("2D SFS computation complete.")

def main():
    parser = argparse.ArgumentParser(
        description="Compute 2D Site Frequency Spectrum (SFS) from VCF.gz files for east and west populations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True, help="Path to input VCF.gz file")
    parser.add_argument("--output", type=str, required=True, help="Path to output 2D SFS file")
    args = parser.parse_args()

    compute_2d_sfs(
        vcf_path=args.vcf,
        output_path=args.output
    )

if __name__ == "__main__":
    main()
