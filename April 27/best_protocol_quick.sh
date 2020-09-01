#!/bin/bash
#SBATCH --job-name=optProt_q

srun julia best_protocol_quick.jl
