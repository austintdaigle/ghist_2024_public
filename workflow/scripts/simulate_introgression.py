import argparse
import msprime
import numpy as np
import gzip
import csv
import demes
import demesdraw
import matplotlib.pyplot as plt

def simulate_neanderthal_admixture(
    Nanc, N_modern, N_neanderthal, N_pop1, N_pop2,
    T_split, T_pop_split, T_admix1, T_admix2,
    m1, m2, mutation_rate, recombination_rate,
    sequence_length=2.5e8
):
    """
    Simulates admixture between two human populations with Neanderthals,
    including variable population sizes and admixture proportions.

    Parameters:
        Nanc (float): Ancestral population size.
        N_modern (float): Effective population size of modern humans.
        N_neanderthal (float): Effective population size of Neanderthals.
        N_pop1 (float): Effective population size of Population 1 after split.
        N_pop2 (float): Effective population size of Population 2 after split.
        T_split (float): Time of split between ancestral and modern humans/Neanderthals.
        T_pop_split (float): Time of split between pop1 and pop2.
        T_admix1 (float): Time of admixture event between Neanderthals and Population 1.
        T_admix2 (float): Time of admixture event between Neanderthals and Population 2.
        m1 (float): Proportion of Neanderthal ancestry in Population 1.
        m2 (float): Proportion of Neanderthal ancestry in Population 2.
        mutation_rate (float): Mutation rate per base per generation.
        recombination_rate (float): Recombination rate per base per generation.
        sequence_length (float): Length of the simulated sequence.

    Returns:
        mutated_ts (tskit.TreeSequence): Tree sequence with simulated mutations.
    """
    print("Simulating Neanderthal admixture with parameters:")
    print(f"  Nanc: {Nanc}, N_modern: {N_modern}, N_neanderthal: {N_neanderthal}")
    print(f"  N_pop1: {N_pop1}, N_pop2: {N_pop2}")
    print(f"  T_split: {T_split}, T_pop_split: {T_pop_split}")
    print(f"  T_admix1: {T_admix1}, T_admix2: {T_admix2}")
    print(f"  m1: {m1}, m2: {m2}")
    print(f"  Mutation rate: {mutation_rate}, Recombination rate: {recombination_rate}")

    # Define demography
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=Nanc)
    demography.add_population(name="modern_humans", initial_size=N_modern)
    demography.add_population(name="neanderthal", initial_size=N_neanderthal)
    
    # Population split events
    # Split "modern_humans" and "neanderthal" from the ancestral population
    demography.add_population_split(
        time=T_split,
        derived=["modern_humans", "neanderthal"],
        ancestral="ancestral"
    )
    # Split "pop1" and "pop2" from "modern_humans" and set initial sizes for each
    demography.add_population(name="pop1", initial_size=N_pop1)
    demography.add_population(name="pop2", initial_size=N_pop2)
    demography.add_population_split(
        time=T_pop_split,
        derived=["pop1", "pop2"],
        ancestral="modern_humans"
    )

    # Admixture events from Neanderthals into Pop1 and Pop2
    demography.add_mass_migration(time=T_admix1, source="pop1", dest="neanderthal", proportion=m1)
    demography.add_mass_migration(time=T_admix2, source="pop2", dest="neanderthal", proportion=m2)

    # Sort events by time to ensure chronological order
    demography.sort_events()

    #try:
    #    # Export demography to demes format for visualization and verification
    #    demes_graph = demography.to_demes()
    #    print("Generating and saving demography plot...")
    #    demesdraw.tubes(demes_graph)
    #    plt.savefig("/nas/longleaf/home/adaigle/ghist_2024/workflow/scripts/demography_plot.png", format="png")
    #    plt.close()
    #except Exception as e:
    #    print(f"Error in converting demography to demes: {e}")
    #    # Optionally, remove or comment out the plotting if demography is invalid
    #    pass

    # Define ancient Neanderthal samples at different times and counts
    samples = [
        msprime.SampleSet(20, population="pop1", time=0),
        msprime.SampleSet(16, population="pop2", time=0),
        msprime.SampleSet(3, population="neanderthal", time=4000),  # AncSite1_1, AncSite1_2, AncSite1_3
        msprime.SampleSet(2, population="neanderthal", time=5000),  # AncSite2_1, AncSite2_2
        msprime.SampleSet(2, population="neanderthal", time=1000),  # AncSite3_1, AncSite3_2
        msprime.SampleSet(2, population="neanderthal", time=1000),  # AncSite4_1, AncSite4_2
        msprime.SampleSet(1, population="neanderthal", time=700)    # AncSite5_1
    ]

    # Simulate ancestry
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=42
    )

    # Simulate mutations
    mutated_ts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        model=msprime.BinaryMutationModel()
    )
    return mutated_ts

def log_simulation_parameters(
    log_file, repnum, Nanc, N_modern, N_neanderthal,
    N_pop1, N_pop2, T_split, T_pop_split,
    T_admix1, T_admix2, m1, m2
):
    """
    Logs the simulation parameters to a specified log file.

    Parameters:
        log_file (str): Path to the log file.
        repnum (int): Replicate number.
        Other parameters: Various simulation parameters to log.
    """
    with open(log_file, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow([
            repnum, Nanc, N_modern, N_neanderthal,
            N_pop1, N_pop2, T_split, T_pop_split,
            T_admix1, T_admix2, m1, m2
        ])
    print(f"Logged parameters for replicate {repnum} to {log_file}")

def main():
    parser = argparse.ArgumentParser(description="Simulate Neanderthal admixture into two human populations.")
    parser.add_argument("--repnum", type=int, required=True, help="Replicate number")
    parser.add_argument("--path", type=str, required=True, help="Path to output file")
    parser.add_argument("--logfile", type=str, default="neanderthal_admixture_log.tsv", help="Path to log file for simulation parameters")
    args = parser.parse_args()

    # Set random parameters for each replicate
    Nanc = np.random.uniform(5000, 50000)
    N_modern = np.random.uniform(5000, 50000)
    N_neanderthal = np.random.uniform(5000, 50000)
    N_pop1 = np.random.uniform(5000, 50000)
    N_pop2 = np.random.uniform(5000, 50000)

    # Ensure T_split > T_pop_split > T_admix1, T_admix2
    T_split = np.random.uniform(10000, 40000)      # e.g., 550,000 years ago in generations
    T_pop_split = np.random.uniform(1000, 3000)    # e.g., 25,000 years ago in generations

    # Ensure T_admix1 and T_admix2 are less than T_pop_split
    if T_pop_split > 1:
        T_admix1 = np.random.uniform(100, T_pop_split - 1)
        T_admix2 = np.random.uniform(100, T_pop_split - 1)
    else:
        T_admix1 = 500
        T_admix2 = 500

    m1 = np.random.uniform(0.01, 0.03)             # Proportion of Neanderthal admixture into Pop1
    m2 = np.random.uniform(0.01, 0.03)             # Proportion of Neanderthal admixture into Pop2

    mutation_rate = 1.29e-8
    recombination_rate = 1.38e-8
    filename = f"{args.path}replicate_{args.repnum}.vcf.gz"

    print("Generated parameters:")
    print(f"  Nanc: {Nanc}, N_modern: {N_modern}, N_neanderthal: {N_neanderthal}")
    print(f"  N_pop1: {N_pop1}, N_pop2: {N_pop2}")
    print(f"  T_split: {T_split}, T_pop_split: {T_pop_split}")
    print(f"  T_admix1: {T_admix1}, T_admix2: {T_admix2}")
    print(f"  m1: {m1}, m2: {m2}")
    print(f"  Mutation rate: {mutation_rate}, Recombination rate: {recombination_rate}")
    print(f"  Output file: {filename}")

    # Log parameters
    log_simulation_parameters(
        args.logfile,
        args.repnum,
        Nanc, N_modern, N_neanderthal,
        N_pop1, N_pop2, T_split, T_pop_split,
        T_admix1, T_admix2, m1, m2
    )

    # Run simulation
    mutated_ts = simulate_neanderthal_admixture(
        Nanc=Nanc,
        N_modern=N_modern,
        N_neanderthal=N_neanderthal,
        N_pop1=N_pop1,
        N_pop2=N_pop2,
        T_split=T_split,
        T_pop_split=T_pop_split,
        T_admix1=T_admix1,
        T_admix2=T_admix2,
        m1=m1,
        m2=m2,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate
    )

    # Modify individual names for VCF output
    individual_names = (
        [f"pop1_{i+1}" for i in range(20)] +
        [f"pop2_{i+1}" for i in range(16)] +
        [f"AncSite1_{i+1}" for i in range(3)] +
        [f"AncSite2_{i+1}" for i in range(2)] +
        [f"AncSite3_{i+1}" for i in range(2)] +
        [f"AncSite4_{i+1}" for i in range(2)] +
        [f"AncSite5_{i+1}" for i in range(1)]
    )

    # Ensure the list length matches the number of individuals in `mutated_ts`
    print(mutated_ts.num_individuals)
    print(len(individual_names))
    if len(individual_names) != mutated_ts.num_individuals:
        raise ValueError("Mismatch between individual_names and the number of individuals in the tree sequence")

    # Write the mutated tree sequence to a gzipped VCF file
    print("Saving to VCF file...")
    with gzip.open(filename, 'wt') as vcf_file:
        mutated_ts.write_vcf(
            vcf_file,
            individual_names=individual_names,
            position_transform=lambda x: [1 + pos for pos in x]  # Shift positions by 1 to make them 1-based
        )

    print(f"VCF file saved successfully to {filename}")

if __name__ == "__main__":
    print("Starting simulation script...")
    main()
    print("Simulation script complete.")
