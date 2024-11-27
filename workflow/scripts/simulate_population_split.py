import argparse
import msprime
import numpy as np
import gzip
import csv

def simulate_population_split(Nanc, N_west, N_east, T, mutation_rate, recombination_rate, sequence_length=1e8):
    """
    Simulates a population split and isolation model with east and west populations.
    
    Parameters:
        Nanc (float): Ancestral population size.
        N_west (float): West population size after split.
        N_east (float): East population size after split (smaller population).
        T (float): Time of population split in generations.
        mutation_rate (float): Mutation rate per base per generation.
        recombination_rate (float): Recombination rate per base per generation.
        sequence_length (float): Length of the simulated sequence.
    
    Returns:
        mutated_ts (tskit.TreeSequence): Mutated tree sequence with population split.
    """
    print("Starting population split simulation with parameters:")
    print(f"  Nanc: {Nanc}, N_west: {N_west}, N_east: {N_east}, T: {T}")
    print(f"  Mutation rate: {mutation_rate}, Recombination rate: {recombination_rate}")
    
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=Nanc)
    demography.add_population(name="west", initial_size=N_west)
    demography.add_population(name="east", initial_size=N_east)

    # Define the population split at time T
    demography.add_population_split(time=T, derived=["east", "west"], ancestral="ancestral")
    print("Demography setup complete. Starting ancestry simulation...")

    # Simulate ancestry with 22 east and 18 west individuals
    ts = msprime.sim_ancestry(
        samples={"west": 18, "east": 22},
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        demography=demography
    )
    print("Ancestry simulation complete. Starting mutation simulation...")

    # Simulate mutations
    mutated_ts = msprime.sim_mutations(ts, rate=mutation_rate, model=msprime.BinaryMutationModel())
    print("Mutation simulation complete.")

    return mutated_ts

def log_simulation_parameters(log_file, repnum, Nanc, N_west, N_east, T):
    """
    Logs the simulation parameters to a specified log file.

    Parameters:
        log_file (str): Path to the log file.
        repnum (int): Replicate number.
        Nanc (float): Ancestral population size.
        N_west (float): West population size after split.
        N_east (float): East population size after split.
        T (float): Time of population split in generations.
    """
    with open(log_file, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow([repnum, Nanc, N_west, N_east, T])
    print(f"Logged parameters for replicate {repnum} to {log_file}")

def main():
    parser = argparse.ArgumentParser(description="Simulate a population split and isolation model.")
    parser.add_argument("--repnum", type=int, required=True, help="Replicate number")
    parser.add_argument("--path", type=str, required=True, help="Path to output file")
    parser.add_argument("--logfile", type=str, default="split_simulation_log.tsv", help="Path to log file for simulation parameters")
    args = parser.parse_args()

    # Set parameters
    Nanc = np.random.uniform(90000, 110000)
    #N_east = np.random.uniform(110000, 130000)
    #N_west = np.random.uniform(18000, 40000)
    N_east = np.random.uniform(15000, 150000)
    N_west = np.random.uniform(15000, 150000)
    T = np.random.uniform(10000, 40000)
    mutation_rate = 5.7e-9
    recombination_rate = 3.386e-9
    filename = f"{args.path}replicate_{args.repnum}.vcf.gz"

    print("Parameters generated:")
    print(f"  Nanc: {Nanc}, N_west: {N_west}, N_east: {N_east}, T: {T}")
    print(f"  Mutation rate: {mutation_rate}, Recombination rate: {recombination_rate}")
    print(f"  Output file: {filename}")

    # Log the parameters
    log_simulation_parameters(
        args.logfile,
        args.repnum,
        Nanc,
        N_west,
        N_east,
        T
    )

    # Run simulation
    mutated_ts = simulate_population_split(
        Nanc=Nanc,
        N_west=N_west,
        N_east=N_east,
        T=T,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate
    )

    # Write the mutated tree sequence to a gzipped VCF file
    print("Saving to VCF file...")
    with gzip.open(filename, 'wt') as vcf_file:
        mutated_ts.write_vcf(
            vcf_file,
            individual_names=[f"east_{i+1}" for i in range(22)] + [f"west_{i+1}" for i in range(18)],
            position_transform=lambda x: [1 + pos for pos in x]  # Shift positions by 1 to make them 1-based
        )

    print(f"VCF file saved successfully to {filename}")

if __name__ == "__main__":
    print("Starting simulation script...")
    main()
    print("Simulation script complete.")