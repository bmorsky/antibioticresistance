#!/usr/bin/env julia

### Load libraries
#using Statistics, OrdinaryDiffEq, DiffEqJump, DiffEqBase, DifferentialEquations, Catalyst
#using StochasticDiffEq, DifferentialEquations, OrdinaryDiffEq, RecursiveArrayTools, LinearAlgebra, DiffEqJump, DiffEqBase

# Test whether this is run on the cluster, or locally
global cluster = any(collect(keys(ENV)).=="SLURM_ARRAY_JOB_ID");

# IF CLUSTER (Used when running on a cluster)
if cluster
	println("You are using SLRUM, your SLURM job ID is ", ENV["SLURM_JOB_ID"])
  global paths="/home/bmorsky/antibiotic"
  global queue = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]);
  # Initialise Parameters
  include(string(paths, "/params.jl"));
  # Initialise Methods
  include(string(paths, "/avgbact.jl")); # in this case LV_dynThresh_spatial.jl
else # (Used when running locally)
  # Initialise Methods
  # global paths="/Users/brycemorsky/"
  include(string(paths, "/methods.jl"));
  # Run Simulation
  for q in 1:25
    global queue = q
	# Initialise Parameters
    include(string(paths, "/params.jl"));
	# Run Simulation
	@time runit();
	println(" with ", queue,".")
  end
end
