# SchriderLab GHIST 2024 Challenge Workflow

This repository contains the pipeline and tools developed to complete the [Genomic History Inference Strategies Tournament (GHIST) 2024](https://gutengroup.arizona.edu/ghist). GHIST is an annual competition that provides simulated population genomic datasets, challenging participants to infer various aspects of the processes that generated those data. The 2024 edition focuses on demographic history inference, offering four challenges of escalating difficulty. We completed the first three, training an MLP to predict the demographic parameters using summary statistics from simulated data. 
<img src="https://github.com/user-attachments/assets/ef6a13c4-ee26-4f8e-9507-bfacfff01fe0" alt="image" width="600"/>

## Overview

We completed some initial parameter estimates in dadi, simulations using msprime, and final estimates using an MLP. 
Relevant scripts for each stage of the analysis are shown below:

### 1. Demographic Inference Using dadi

- **Objective**: Establish priors through basic demographic models. 
- **Implementation**: Jupyter notebooks for specific models are located in `workflow/scripts/dadi_<model_name>.ipynb`.

### 2. Reproducible Workflow with Snakemake

- **Objective**: Perform population genetics simulations and save outputs in VCF format, then compute summary statistics using pylibseq and scikit-allel.
- **Implementation**: The Snakemake pipeline is defined in `workflow/Snakefile`, which will need to be modified to specify which scenario you want to simulate. Scripts for various scenarios (e.g., bottleneck, population split) are available in `workflow/scripts`.

#### Running the Pipeline

To execute the pipeline, use the following command:

```bash
snakemake -d <output_dir>
```
Replace <output_dir> with the path to your output directory.

### 3. Machine Learning Pipeline

- **Objective**: Train neural networks and conduct experiments for demographic model classification.
- **Implementation**: Jupyter notebooks for MLP training along with some experiments and analyses are found in `workflow/scripts/mlp_<model>.ipynb`.

## Environment Setup

To ensure compatibility and reproducibility, set up the Conda environment as follows:

### Step 1: Create the Environment

```bash
conda env create -f workflow/envs/ghist.yaml
```

### Step 2: Activate the Environment

```bash
conda activate ghist
```
