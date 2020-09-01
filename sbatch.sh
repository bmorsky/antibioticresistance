#!/bin/bash
#SBATCH --job-name=antibiotic
#SBATCH --array=1-50%50

srun julia runSims.jl
