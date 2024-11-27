import argparse
import numpy as np
import pandas as pd
import allel
import sys

def compute_2d_sfs(vcf_path, output_path):
    print(f"Loading VCF file from {vcf_path}...")

    try:
        # Load genotype data and sample names from VCF
        callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'samples'])
        genotypes = allel.GenotypeArray(callset['calldata/GT'])  # Shape: (variants, samples, ploidy)
        sample_names = callset['samples']
    
        # Separate indices for island and mainland samples
        island_indices = [i for i, name in enumerate(sample_names) if name.startswith('island')]
        mainland_indices = [i for i, name in enumerate(sample_names) if name.startswith('mainland')]
    
        # Genotypes for each population
        genotypes_island = genotypes[:, island_indices]
        genotypes_mainland = genotypes[:, mainland_indices]
    
    except Exception as e:
        print(f"Error reading VCF file: {e}")
        sys.exit(1)

    print("Computing 2D Site Frequency Spectrum (SFS) for island and mainland populations...")

    # Compute derived allele counts for each population
    ac_island = genotypes_island.count_alleles()
    ac_mainland = genotypes_mainland.count_alleles()
    dac_island = ac_island[:, 1]  # Derived allele counts for island
    dac_mainland = ac_mainland[:, 1]  # Derived allele counts for mainland

    # Compute 2D SFS
    sfs_2d = allel.joint_sfs(dac_island, dac_mainland)

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
        description="Compute 2D Site Frequency Spectrum (SFS) from VCF.gz files for island and mainland populations.",
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
