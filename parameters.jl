### 1- Initialize parameters ...
#println("Here is your queue number: ", queue)
grid = collect(Base.product(
# lines with value = 99 are currently depracted
[1e-4 1e-5 1e-6 1e-7], #1 inverse carrying capacity
[0.1 1 10 100], #2 competition parameter κ
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
global γ = grid[1] # inverse carrying capacity
global κ = grid[2] # competition parameter κ

global num_sims = 2 # number of realisations
global time = 100 # length of a simulation
global X₀ = 1e4 # initial number of susceptible, X

# Parameters for step lengths and periods
global step_length = 5
global num_steps = 5
global period = step_length*num_steps
