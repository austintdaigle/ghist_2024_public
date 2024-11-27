import argparse
import numpy as np
import pandas as pd
import allel
import scipy.stats
import sys
import libsequence

def compute_summary_stats(vcf_path, window_length, output_path, max_dist_ld=1000, bin_size_ld=100, sequence_length=1e8):
    print(f"Loading VCF file from {vcf_path}...")

    try:
        # Load genotype data, variant positions, and sample names from VCF
        callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'variants/POS', 'samples'])
        genotypes = allel.GenotypeArray(callset['calldata/GT'])  # Shape: (variants, samples, ploidy)
        positions = np.array(callset['variants/POS'], dtype=float)
        sample_names = callset['samples']
    
        # Separate indices for island and mainland samples
        island_indices = [i for i, name in enumerate(sample_names) if name.startswith('island')]
        mainland_indices = [i for i, name in enumerate(sample_names) if name.startswith('mainland')]
    
        # Genotypes for each population
        genotypes_joint = genotypes  # For clarity
        genotypes_island = genotypes[:, island_indices]
        genotypes_mainland = genotypes[:, mainland_indices]
    
    except Exception as e:
        print(f"Error reading VCF file: {e}")
        sys.exit(1)

    # Compute summary statistics for joint, island, and mainland populations
    stats_joint = compute_population_stats(
        genotypes=genotypes_joint,
        positions=positions,
        label='joint',
        max_dist_ld=max_dist_ld,
        bin_size_ld=bin_size_ld,
        window_length=window_length,
        sequence_length=sequence_length
    )

    stats_island = compute_population_stats(
        genotypes=genotypes_island,
        positions=positions,
        label='island',
        max_dist_ld=max_dist_ld,
        bin_size_ld=bin_size_ld,
        window_length=window_length,
        sequence_length=sequence_length
    )

    stats_mainland = compute_population_stats(
        genotypes=genotypes_mainland,
        positions=positions,
        label='mainland',
        max_dist_ld=max_dist_ld,
        bin_size_ld=bin_size_ld,
        window_length=window_length,
        sequence_length=sequence_length
    )

    # Compute pairwise statistics between island and mainland populations
    stats_pairwise = compute_pairwise_stats(
        genotypes1=genotypes_island,
        genotypes2=genotypes_mainland,
        positions=positions,
        label='island_mainland'
    )

    # Combine all statistics into a single dictionary
    all_stats = {}
    all_stats.update(stats_joint)
    all_stats.update(stats_island)
    all_stats.update(stats_mainland)
    all_stats.update(stats_pairwise)

    df = pd.DataFrame([all_stats])
    print(f"Saving summary statistics to {output_path}...")
    try:
        df.to_csv(output_path, index=False, sep='\t')
    except Exception as e:
        print(f"Error saving summary statistics: {e}")
        sys.exit(1)

    print("Summary statistics computation complete.")

def compute_population_stats(genotypes, positions, label, max_dist_ld, bin_size_ld, window_length, sequence_length):
    """
    Computes summary statistics for a given population.
    """
    print(f"Computing summary statistics for {label} population...")

    # Convert genotypes to allele counts
    gn = genotypes.to_n_alt()  # Shape: (variants, samples)

    # Compute LD decay
    genotype_counts = gn
    ld_decay_stats = compute_ld_decay(genotype_counts, positions, max_dist=max_dist_ld, bin_size=bin_size_ld)

    # Prepare data for libsequence
    haplotypes = genotypes.to_haplotypes().tolist()
    positions_list = positions.tolist()

    # Create VariantMatrix
    vm = libsequence.VariantMatrix(haplotypes, positions_list)

    # Compute Runs of Homozygosity (ROH)
    print("Computing Runs of Homozygosity (ROH)...")
    # Define minimum ROH length in base pairs
    min_roh_length = 10000  # 10kb

    # Convert positions to integer array
    pos_int = positions.astype(int)
    # Calculate genotype mask: True if homozygous (0 or 2), False otherwise
    homozygous = (gn == 0) | (gn == 2)

    # Initialize list to store fROH for each individual
    fROH_list = []

    for individual_homo in homozygous.T:
        # Identify runs of homozygosity
        is_roh = individual_homo.astype(int)
        roh_starts = np.where(np.diff(np.concatenate(([0], is_roh))) == 1)[0]
        roh_ends = np.where(np.diff(np.concatenate((is_roh, [0]))) == -1)[0]

        # Calculate ROH lengths
        roh_lengths = []
        for start, end in zip(roh_starts, roh_ends):
            length = pos_int[end - 1] - pos_int[start]
            if length >= min_roh_length:
                roh_lengths.append(length)

        # Calculate fROH for the individual
        total_roh_length = np.sum(roh_lengths)
        fROH = total_roh_length / sequence_length
        fROH_list.append(fROH)

    # Calculate mean, median, SD, and IQR of fROH
    mean_fROH = np.mean(fROH_list)
    median_fROH = np.median(fROH_list)
    std_fROH = np.std(fROH_list)
    iqr_fROH = scipy.stats.iqr(fROH_list)

    print(f"Mean fROH for {label}: {mean_fROH}")

    # Existing statistics computations
    print("Computing summary statistics in windows...")
    bin_size = window_length  # Use provided window length
    thetapi_values = []
    hprime_values = []
    thetaw_values = []
    faywuh_values = []
    tajd_values = []
    haplotype_diversity_values = []
    number_of_haplotypes_values = []

    # Calculate statistics for each window
    lefts = range(0, int(sequence_length), bin_size)
    for l in lefts:
        # Get windowed VariantMatrix
        w = vm.window(l, l + bin_size)
        if w.nsites == 0:
            continue  # Skip windows with no variants

        acw = w.count_alleles()
        thetapi_value = libsequence.thetapi(acw) / bin_size
        hprime_value = libsequence.hprime(acw, 0)
        thetaw_value = libsequence.thetaw(acw) / bin_size
        faywuh_value = libsequence.faywuh(acw, 0)
        tajd_value = libsequence.tajd(acw)

        hap_div_value = libsequence.haplotype_diversity(w)
        num_hap_value = libsequence.number_of_haplotypes(w)

        if not np.isnan(thetapi_value) and not np.isnan(hprime_value):
            thetapi_values.append(thetapi_value)
            hprime_values.append(hprime_value)
            thetaw_values.append(thetaw_value)
            faywuh_values.append(faywuh_value)
            tajd_values.append(tajd_value)
            haplotype_diversity_values.append(hap_div_value)
            number_of_haplotypes_values.append(num_hap_value)

    # Compile statistics
    stats = {
        f'{label}_thetapi_mean': np.mean(thetapi_values),
        f'{label}_thetapi_median': np.median(thetapi_values),
        f'{label}_thetapi_sd': np.std(thetapi_values),
        f'{label}_thetapi_iqr': scipy.stats.iqr(thetapi_values),
        f'{label}_hprime_mean': np.mean(hprime_values),
        f'{label}_hprime_median': np.median(hprime_values),
        f'{label}_hprime_sd': np.std(hprime_values),
        f'{label}_hprime_iqr': scipy.stats.iqr(hprime_values),
        f'{label}_thetaw_mean': np.mean(thetaw_values),
        f'{label}_thetaw_median': np.median(thetaw_values),
        f'{label}_thetaw_sd': np.std(thetaw_values),
        f'{label}_thetaw_iqr': scipy.stats.iqr(thetaw_values),
        f'{label}_faywuh_mean': np.mean(faywuh_values),
        f'{label}_faywuh_median': np.median(faywuh_values),
        f'{label}_faywuh_sd': np.std(faywuh_values),
        f'{label}_faywuh_iqr': scipy.stats.iqr(faywuh_values),
        f'{label}_tajd_mean': np.mean(tajd_values),
        f'{label}_tajd_median': np.median(tajd_values),
        f'{label}_tajd_sd': np.std(tajd_values),
        f'{label}_tajd_iqr': scipy.stats.iqr(tajd_values),
        f'{label}_haplotype_diversity_mean': np.mean(haplotype_diversity_values),
        f'{label}_haplotype_diversity_median': np.median(haplotype_diversity_values),
        f'{label}_haplotype_diversity_sd': np.std(haplotype_diversity_values),
        f'{label}_haplotype_diversity_iqr': scipy.stats.iqr(haplotype_diversity_values),
        f'{label}_number_of_haplotypes_mean': np.mean(number_of_haplotypes_values),
        f'{label}_number_of_haplotypes_median': np.median(number_of_haplotypes_values),
        f'{label}_number_of_haplotypes_sd': np.std(number_of_haplotypes_values),
        f'{label}_number_of_haplotypes_iqr': scipy.stats.iqr(number_of_haplotypes_values),
        f'{label}_fROH_mean': mean_fROH,
        f'{label}_fROH_median': median_fROH,
        f'{label}_fROH_sd': std_fROH,
        f'{label}_fROH_iqr': iqr_fROH
    }

    print("Computing Site Frequency Spectrum (SFS)...")
    # Compute allele counts per variant
    ac = genotypes.count_alleles()

    # Compute derived allele counts (assuming allele '1' is derived)
    dac = ac[:, 1]  # 1D array of derived allele counts per variant

    # Set the total number of chromosomes (alleles) for diploid data
    n = genotypes.shape[1] * 2  # Total number of chromosomes

    # Compute SFS
    sfs = allel.sfs(dac, n=n)

    # Add SFS to stats
    for i, sfs_col in enumerate(sfs):
        stats[f'{label}_sfs_{i+1}'] = sfs_col

    # Adding LD decay stats, prefixing keys with label
    ld_decay_stats_prefixed = {f'{label}_{k}': v for k, v in ld_decay_stats.items()}
    stats.update(ld_decay_stats_prefixed)

    # Add sample and variant counts
    stats[f'{label}_num_samples'] = gn.shape[1]
    stats[f'{label}_num_variants'] = gn.shape[0]

    return stats

def compute_pairwise_stats(genotypes1, genotypes2, positions, label):
    """
    Computes statistics comparing two populations.
    """
    print(f"Computing pairwise statistics between {label} populations...")

    # Compute allele counts per variant
    ac1 = genotypes1.count_alleles()
    ac2 = genotypes2.count_alleles()

    # Compute Hudson's Fst
    num, den = allel.hudson_fst(ac1, ac2)
    # Handle division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        fst_per_site = np.divide(num, den)
    # Remove NaNs and Infs
    fst_per_site = fst_per_site[np.isfinite(fst_per_site)]
    # Compute mean Fst
    if np.sum(den) > 0:
        mean_fst = np.sum(num) / np.sum(den)
    else:
        mean_fst = np.nan

    # Compute Weir and Cockerham's Fst
    genotypes_combined = allel.GenotypeArray(
        np.concatenate([genotypes1, genotypes2], axis=1))
    subpops = [
        list(range(genotypes1.n_samples)),
        list(range(genotypes1.n_samples, genotypes_combined.n_samples))
    ]

    # Compute components for Weir and Cockerham's Fst
    a, b, c = allel.weir_cockerham_fst(genotypes_combined, subpops)
    # Compute per-site Fst
    with np.errstate(divide='ignore', invalid='ignore'):
        fst_wc_per_site = np.divide(a, a + b + c)
    fst_wc_per_site = fst_wc_per_site[np.isfinite(fst_wc_per_site)]
    # Compute mean Weir and Cockerham's Fst
    sum_a = np.sum(a)
    sum_b = np.sum(b)
    sum_c = np.sum(c)
    if sum_a + sum_b + sum_c > 0:
        mean_fst_wc = sum_a / (sum_a + sum_b + sum_c)
    else:
        mean_fst_wc = np.nan

    # Compute Patterson's F2 statistic
    f2 = allel.patterson_f2(ac1, ac2)
    mean_f2 = np.nanmean(f2)

    stats = {
        f'{label}_fst_hudson_mean': mean_fst,
        f'{label}_fst_weir_cockerham_mean': mean_fst_wc,
        f'{label}_f2_mean': mean_f2
    }

    return stats

def compute_ld_decay(gn, positions, max_dist=1000, bin_size=100):
    """
    Computes LD decay using Rogers and Huff's r^2 within a maximum distance.
    """
    print(f"Computing LD decay within {max_dist} bp and bin size {bin_size} bp...")
    num_bins = int(np.ceil(max_dist / bin_size))
    ld_decay = np.zeros(num_bins)
    ld_counts = np.zeros(num_bins)
    bin_edges = np.arange(0, max_dist + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    num_variants = gn.shape[0]
    for i in range(num_variants):
        start_pos = positions[i]
        # Find the index where positions exceed start_pos + max_dist
        max_index = np.searchsorted(positions, start_pos + max_dist, side='right')
        # Ensure we have at least two SNPs to compute LD
        if max_index > i + 1:
            # Extract genotype window for current SNP and subsequent SNPs within max_dist
            gn_window = gn[i:max_index]
            # Compute Rogers and Huff's r
            try:
                r_values = allel.rogers_huff_r(gn_window)
                r2_values = np.square(r_values)
                # Ensure r2_values is an array, even if only one value is present
                if np.isscalar(r2_values):
                    r2_values = np.array([r2_values])
            except Exception as e:
                continue

            # Compute pairwise distances
            pos_window = positions[i:max_index]
            dists = pos_window - start_pos  # since pos_window starts at start_pos

            # Bin r^2 values by distance
            bin_indices = np.digitize(dists, bin_edges) - 1

            # Iterate over r^2 and bin_idx
            for bin_idx, r2 in zip(bin_indices, r2_values):
                if 0 <= bin_idx < num_bins:
                    if not np.isnan(r2) and not np.isinf(r2):
                        ld_decay[bin_idx] += r2
                        ld_counts[bin_idx] += 1

    if np.all(ld_counts == 0):
        print("No SNP pairs found within the specified distance for LD decay calculation.")
        ld_decay_stats = {f'ld_decay_{int(center)}bp': 0 for center in bin_centers}
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            ld_decay_mean = np.divide(ld_decay, ld_counts)
            ld_decay_mean = np.nan_to_num(ld_decay_mean)  # Replace NaNs with 0

        ld_decay_stats = {f'ld_decay_{int(center)}bp': mean_r2 for center, mean_r2 in zip(bin_centers, ld_decay_mean)}

    return ld_decay_stats

def main():
    parser = argparse.ArgumentParser(
        description="Compute summary statistics from VCF.gz files with specified window length.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True, help="Path to input VCF.gz file")
    parser.add_argument("--window_length", type=int, required=True, help="Window length in base pairs for summary statistics")
    parser.add_argument("--output", type=str, required=True, help="Path to output summary statistics file")
    parser.add_argument("--max_dist_ld", type=int, default=1000, help="Maximum distance in base pairs for LD decay")
    parser.add_argument("--bin_size_ld", type=int, default=100, help="Bin size in base pairs for LD decay")
    parser.add_argument("--sequence_length", type=float, default=1e8, help="Total sequence length")
    args = parser.parse_args()

    compute_summary_stats(
        vcf_path=args.vcf,
        window_length=args.window_length,
        output_path=args.output,
        max_dist_ld=args.max_dist_ld,
        bin_size_ld=args.bin_size_ld,
        sequence_length=args.sequence_length
    )

if __name__ == "__main__":
    main()
