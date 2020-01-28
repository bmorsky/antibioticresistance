#!/bin/bash
#SBATCH --job-name=tpbvp

srun julia tpbvp.jl
