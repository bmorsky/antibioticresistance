#!/bin/bash
#SBATCH --job-name=antibiotic
#SBATCH --array=1-32%6

srun julia alpha.jl
srun julia kappa.jl
srun julia mu.jl
srun julia gamma.jl
srun julia xi.jl
