#!/bin/bash
#SBATCH --job-name=best_protocol

srun julia gillespie_avg_kim_2020.jl
