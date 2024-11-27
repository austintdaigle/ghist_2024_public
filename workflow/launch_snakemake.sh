#!/bin/bash
unset SLURM_MEM_PER_CPU
unset SLURM_MEM_PER_NODE
unset SLURM_MEM_PER_GPU

# Now call Snakemake
snakemake -d /nas/longleaf/home/adaigle/work/ghist_2024_work/bottleneck --executor slurm --workflow-profile /nas/longleaf/home/adaigle/ghist_2024/workflow/profiles/slurm/