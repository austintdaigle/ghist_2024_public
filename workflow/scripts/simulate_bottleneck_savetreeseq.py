import os
import argparse
import msprime
import tskit
import random
import gzip
import csv

def simulate_bottleneck(initial_population_size, bottleneck_population_size, bottleneck_time, sequence_length, recombination_rate, mutation_rate):
    """
    Simulates a bottleneck event and returns the mutated tree sequence.

    Parameters:
        initial_population_size (int): Initial population size.
        bottleneck_population_size (int): Population size during the bottleneck.
        bottleneck_time (int): Time (in generations) when the bottleneck occurs.
        sequence_length (float): Length of the simulated sequence.
        recombination_rate (float): Recombination rate per base per generation.
        mutation_rate (float): Mutation rate per base per generation.

    Returns:
        mutated_ts (tskit.TreeSequence): The mutated tree sequence.
    """
    bottleneck_fraction = bottleneck_population_size / initial_population_size
    print(f"Simulating a bottleneck to this proportion of initial size: {bottleneck_fraction}")

    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=bottleneck_population_size)
    demography.add_population_parameters_change(time=bottleneck_time, initial_size=initial_population_size)

    # Simulate ancestry
    ts = msprime.sim_ancestry(
        samples=20,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        demography=demography,
    )

    # Simulate mutations using a binary mutation model
    mutated_ts = msprime.sim_mutations(ts, rate=mutation_rate, model=msprime.BinaryMutationModel())

    return mutated_ts

def log_simulation_parameters(log_file, repnum, initial_population_size, bottleneck_population_size, bottleneck_time):
    """
    Logs the simulation parameters to a specified log file.

    Parameters:
        log_file (str): Path to the log file.
        repnum (int): Replicate number.
        initial_population_size (int): Initial population size.
        bottleneck_population_size (int): Population size during the bottleneck.
        bottleneck_time (int): Time (in generations) when the bottleneck occurs.
    """
    with open(log_file, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow([repnum, initial_population_size, bottleneck_population_size, bottleneck_time])
    print(f"Logged parameters for replicate {repnum} to {log_file}")

def main():
    parser = argparse.ArgumentParser(description="Simulate a population bottleneck and save as VCF.")
    parser.add_argument("--repnum", type=int, required=True, help="Replicate number")
    parser.add_argument("--path", type=str, required=True, help="Path to output file")
    parser.add_argument("--logfile", type=str, default="simulation_log.tsv", help="Path to log file for simulation parameters")
    args = parser.parse_args()

    # Randomly generate population parameters
    initial_population_size = random.randint(10_000, 15_000)
    bottleneck_population_size = random.randint(10, 2_000)
    bottleneck_time = random.randint(10, 1_000)

    filename = f"{args.path}replicate_{args.repnum}.vcf.gz"

    # Hardcoded values for sequence length, recombination rate, and mutation rate
    sequence_length = 1e8
    recombination_rate = 1.007e-8
    mutation_rate = 1.26e-8

    # Simulate the bottleneck and get the tree sequence
    mutated_ts = simulate_bottleneck(
        initial_population_size,
        bottleneck_population_size,
        bottleneck_time,
        sequence_length,
        recombination_rate,
        mutation_rate
    )

    # Write the mutated tree sequence directly to a gzipped VCF file
    with gzip.open(filename, 'wt') as vcf_file:
        mutated_ts.write_vcf(
            vcf_file,
            position_transform=lambda x: [1 + pos for pos in x]  # Shift positions by 1 to make them 1-based
        )

    print(f"VCF file saved directly to {filename}")

    # Log the parameters to a file
    log_simulation_parameters(
        args.logfile,
        args.repnum,
        initial_population_size,
        bottleneck_population_size,
        bottleneck_time
    )

if __name__ == "__main__":
    main()