#!/bin/bash
#SBATCH --job-name=best_protocol

srun julia best_protocol.jl
