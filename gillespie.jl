# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
using DiffEqBiological, DifferentialEquations, Distributions, Plots
# Parameters
b₁ = 0.115
d₁ = 0.01
b₂ = 0.111
d₂ = 0.01
γ₁₁ = 0.0001
γ₁₂ = 0.00005
γ₂₁ = 0.0002
γ₂₂ = 0.0001
α₁ = 0.187
α₂ = 0.002
μₐ = 0.00015
μ = 0.000015
ζ = 1.0
η = 0.000001

# Initial conditions
control = 2.0 # initial input rate of bacteria
time = 100 # length of process
init_pop = 1000.0 # initial number of susceptible, n₁
rates = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, μₐ, ζ, η, control)

# Protocol parameters
durOn = 2.0 #duration the antibiotic is off
durOff = 2.0 #duration the antibiotic is on
num_periods = 10 #number of periods
num_turns = num_periods*(durOn+durOff) #number of turns

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
    b₁, n₁ → 2n₁ # n₁ birth
    d₁, n₁ → 0 # n₁ death
	b₂, n₂ → 2n₂ # n₂ birth
	d₂, n₂ → 0 # n₂ death
	γ₁₁, n₁ + n₁ → 0 # self competition of n₁
	γ₁₂, n₁ + n₂ → 0 # competition to n₁ from n₂
	γ₂₁, n₂ + n₁ → 0 # competition to n₂ from n₁
	γ₂₂, n₂ + n₂ → 0 # self competition of n₂
	α₁, a + n₁ → a # death by antibiotic
	α₂, a + n₂ → a # death by antibiotic, α₂ < α₁
	μ, n₁ → n₂ # mutation
	μ, n₂ → n₁ # mutation
	μₐ, a + n₁ → a + n₂ # mutation due to stress from the antibiotic
	ζ, a → 0 # dissipation of the antibiotic
	η, n₁ + n₂ → 2n₂ # plasmid (resistance) transfer
	control, 0 → a # inflow of antibiotic
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ ζ η control

# Controller conditions
function condition(t,u,integrator)
	mod(integrator.t,10)*mod(integrator.t,5) # condition when dynamics change
end
function affect!(t,u,integrator) # effect that occurs when the condition is met
	if integrator.u == 2
		integrator.u = 0
	else
		integrator.u = 2
	end
end
cb = ContinuousCallback(condition,affect!)

# Solve system 40 times and plot the results
for sim = 1:40
	prob = DiscreteProblem([init_pop; 0; control],(0.0,time),rates)
	jump_prob = JumpProblem(prob,Direct(),LV_model)
	sol = solve(jump_prob,FunctionMap(),callback=cb)
	plot!(sol)
end
plot!(legend=false)
