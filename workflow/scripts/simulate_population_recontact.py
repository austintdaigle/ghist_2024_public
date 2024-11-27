import argparse
import msprime
import numpy as np
import gzip
import csv
import demes  # Import demes for demography export
import demesdraw  # Import demesdraw for visualization
import matplotlib.pyplot as plt  # Import matplotlib for saving the plot

def simulate_population_split_secondary_contact(Nanc, N_mainland, N_island, T_split, T_contact, m, mutation_rate, recombination_rate, sequence_length=1e8):
    """
    Simulates a population split with initial isolation, followed by secondary contact with one-way migration from mainland to island.

    Parameters:
        Nanc (float): Ancestral population size.
        N_mainland (float): Mainland population size after split.
        N_island (float): Island population size after split.
        T_split (float): Time of initial population split in generations.
        T_contact (float): Time of secondary contact and migration onset.
        m (float): Migration rate from mainland to island.
        mutation_rate (float): Mutation rate per base per generation.
        recombination_rate (float): Recombination rate per base per generation.
        sequence_length (float): Length of the simulated sequence.

    Returns:
        mutated_ts (tskit.TreeSequence): Mutated tree sequence with split and one-way secondary contact migration.
    """
    print("Starting population split with one-way secondary contact migration simulation with parameters:")
    print(f"  Nanc: {Nanc}, N_mainland: {N_mainland}, N_island: {N_island}, T_split: {T_split}, T_contact: {T_contact}, m: {m}")
    print(f"  Mutation rate: {mutation_rate}, Recombination rate: {recombination_rate}")

    # Define demography
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=Nanc)
    demography.add_population(name="mainland", initial_size=N_mainland)
    demography.add_population(name="island", initial_size=N_island)

    # Define the population split at time T_split
    demography.add_population_split(time=T_split, derived=["mainland", "island"], ancestral="ancestral")

    # Set one-way migration from mainland to island with rate m starting at time 0
    demography.add_migration_rate_change(time=0.0, rate=m, source="island", dest="mainland")

    # Stop migration at T_contact by setting the migration rate to 0
    demography.add_migration_rate_change(time=T_contact, rate=0.0, source="island", dest="mainland")

    # Sort events by time to ensure chronological order
    demography.sort_events()

    # Configure samples with diploid ploidy
    samples = [
        msprime.SampleSet(22, population="mainland", ploidy=2, time=0),
        msprime.SampleSet(18, population="island", ploidy=2, time=0)
    ]
       #    # Export demography to demes format for visualization and verification
    #demes_graph = demography.to_demes()
    ## Plot and save the demography as an image
    #print("Generating and saving demography plot...")
    #demesdraw.tubes(demes_graph)
    #plt.savefig("/nas/longleaf/home/adaigle/ghist_2024/workflow/scripts/recontactmodel.png", format="png")
    #plt.close()
    ##print(f"Demography plot saved as {image_output}")
    #print("Demography setup complete. Starting ancestry simulation...")
    #exit()
    # Simulate ancestry with 22 mainland and 8 island individuals
    ts = msprime.sim_ancestry(
        samples={"mainland": 22, "island": 8},
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        demography=demography
    )
    print("Ancestry simulation complete. Starting mutation simulation...")

    # Simulate mutations
    mutated_ts = msprime.sim_mutations(ts, rate=mutation_rate, model=msprime.BinaryMutationModel())
    print("Mutation simulation complete.")

    return mutated_ts

def log_simulation_parameters(log_file, repnum, Nanc, N_mainland, N_island, T_split, T_contact, m):
    """
    Logs the simulation parameters to a specified log file.

    Parameters:
        log_file (str): Path to the log file.
        repnum (int): Replicate number.
        Nanc (float): Ancestral population size.
        N_mainland (float): Mainland population size after split.
        N_island (float): Island population size after split.
        T_split (float): Time of population split in generations.
        T_contact (float): Time of secondary contact and migration onset.
        m (float): Migration rate between island and mainland.
    """
    with open(log_file, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow([repnum, Nanc, N_mainland, N_island, T_split, T_contact, m])
    print(f"Logged parameters for replicate {repnum} to {log_file}")

def main():
    parser = argparse.ArgumentParser(description="Simulate a population split with secondary contact model.")
    parser.add_argument("--repnum", type=int, required=True, help="Replicate number")
    parser.add_argument("--path", type=str, required=True, help="Path to output file")
    parser.add_argument("--logfile", type=str, default="split_secondary_contact_log.tsv", help="Path to log file for simulation parameters")
    args = parser.parse_args()

    # Set parameters based on specified ranges
    Nanc = np.random.uniform(12000, 120000)
    N_mainland = np.random.uniform(25000, 250000)
    N_island = np.random.uniform(5000, 50000)
    T_split = np.random.uniform(10000, 30000)
    T_contact = np.random.uniform(10, 1000)
    m = 10 ** np.random.uniform(-5, -2)
    mutation_rate = 5.7e-9
    recombination_rate = 3.386e-9
    filename = f"{args.path}replicate_{args.repnum}.vcf.gz"

    print("Parameters generated:")
    print(f"  Nanc: {Nanc}, N_mainland: {N_mainland}, N_island: {N_island}")
    print(f"  T_split: {T_split}, T_contact: {T_contact}, m: {m}")
    print(f"  Mutation rate: {mutation_rate}, Recombination rate: {recombination_rate}")
    print(f"  Output file: {filename}")

    # Log the parameters
    log_simulation_parameters(
        args.logfile,
        args.repnum,
        Nanc,
        N_mainland,
        N_island,
        T_split,
        T_contact,
        m
    )

    # Run simulation
    mutated_ts = simulate_population_split_secondary_contact(
        Nanc=Nanc,
        N_mainland=N_mainland,
        N_island=N_island,
        T_split=T_split,
        T_contact=T_contact,
        m=m,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate
    )

    # Write the mutated tree sequence to a gzipped VCF file
    print("Saving to VCF file...")
    with gzip.open(filename, 'wt') as vcf_file:
        mutated_ts.write_vcf(
            vcf_file,
            individual_names=[f"mainland_{i+1}" for i in range(22)] + [f"island_{i+1}" for i in range(8)],
            position_transform=lambda x: [1 + pos for pos in x]  # Shift positions by 1 to make them 1-based
        )

    print(f"VCF file saved successfully to {filename}")

if __name__ == "__main__":
    print("Starting simulation script...")
    main()
    print("Simulation script complete.")
