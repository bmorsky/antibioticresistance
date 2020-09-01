#!/bin/bash
#SBATCH --job-name=gillespie_avg

srun julia gillespie_avg_kim_nov30.jl
