### 1- Initialize parameters ...
#println("Here is your queue number: ", queue)
mygrid = collect(Base.product(
# lines with value = 99 are currently depracted
[4.0], #1 competition parameter κ
[1e-4,1e-5,1e-6], #2 mutation rate μ
[1.0,10.0], #3 mutant stress m = μ′/μ
[0], #4 stochasticity
[1e-5], #5 interaction rate γ
[0.4] #6 antibiotic kill rate α
))[queue]; # additional repetitions (global; to spread jobs over more cores)
# # Record parameters for export
# setup = Dict(
#   # Setup for simulation
#   "nRound"            => grid[1], # ...
#   );
#
# # Make parameters globally accessable
# nRound            = setup["nRound"];

## Set local parameters
global κ = mygrid[1] # competition parameter κ
global μ = mygrid[2] # mutation rate μ
global m = mygrid[3] # mutant stress σ
global ξ = mygrid[4] # stochasticity ξ
global γ = mygrid[5] # interaction rate γ
global α = mygrid[6] # antibiotic kill rate α

global num_sims = 50 # number of realisations
global fin_time = 1344 #336# length of a simulation
#global X₀ = 1e5 # initial number of susceptible, X

# Parameters for step lengths and periods
global step_length = 4 #1 #4
global num_steps = 12 #12
global period = step_length*num_steps
