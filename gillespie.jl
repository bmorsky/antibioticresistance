# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
using DiffEqBiological, DifferentialEquations, Distributions, Plots, ParameterizedFunctions
# Parameters
b₁ = 0.02
d₁ = 0.005
b₂ = 0.02
d₂ = 0.005
γ₁₁ = 1.5*1e-5
γ₁₂ = 1.5*0.88e-5
γ₂₁ = 1.5*1.136e-5
γ₂₂ = 1.5*1e-5
α₁ = 0.0015
α₂ = 0.02*α₁
μ = 1e-6
ζ = 0.0
η = 0 #4e-3
controlrate = 0
# Initial conditions
control = 1000.0 # initial input rate of bacteria
time = 1000 # length of process
init_pop = 1e3 # initial number of susceptible, n₁
rates = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, ζ, η, controlrate)

# # Protocol parameters
# durOn = 2.0 #duration the antibiotic is off
# durOff = 2.0 #duration the antibiotic is on
# num_periods = 10 #number of periods
# num_turns = num_periods*(durOn+durOff) #number of turns

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
	10*μ, a + n₁ → a + n₂ # mutation due to stress from the antibiotic
	ζ, a → 0 # dissipation of the antibiotic
	η, n₁ + n₂ → 2n₂ # plasmid (resistance) transfer
	controlrate, 0 → a # inflow of antibiotic
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ ζ η controlrate

# Controller conditions
protocolOn = collect(0:80:1000)
protocolOff = collect(79:80:1000)
protocol = merge(protocolOn,protocolOff)
protocol = collect(0:1:time)
function conditionOn(u,t,integrator)
	t in (protocolOn)
end
function conditionOff(u,t,integrator)
	t in (protocolOff)
end
function conditionStop(u,t,integrator)
	u[1]+u[2] == 0
end
function affectOn!(integrator) # effect that occurs when the condition is met
	#integrator.u[3] += 200
	integrator.u[3] = control
end
function affectOff!(integrator) # effect that occurs when the condition is met
	integrator.u[3] = 0
end
function affectStop!(integrator) # effect that occurs when the condition is met
	terminate!(integrator)
end
cbOn = DiscreteCallback(conditionOn,affectOn!)
cbOff = DiscreteCallback(conditionOff,affectOff!)
cbStop = DiscreteCallback(conditionStop,affectStop!)
cbs = CallbackSet(cbOn,cbOff,cbStop)
prob = DiscreteProblem([init_pop; 0; control],(0.0,time),rates)
jump_prob = JumpProblem(prob,Direct(),LV_model)
sol = solve(jump_prob,FunctionMap(),callback=cbs,tstops=protocol)
plot(sol)

# # Solve system 40 times and plot the results
# for sim = 1:1
# 	prob = DiscreteProblem([init_pop; 0; 10],(0.0,time),rates)
# 	jump_prob = JumpProblem(prob,Direct(),LV_model)
# 	sol = solve(jump_prob,FunctionMap(),callback=cb)
# 	plot!(sol)
# end
# plot!(legend=false)









# # Solve system 40 times and plot the results
# for sim = 1:1
# 	prob = DiscreteProblem([init_pop; 0; 10],(0.0,time),rates)
# 	jump_prob = JumpProblem(prob,Direct(),LV_model)
# 	sol = solve(jump_prob,FunctionMap(),callback=cb)
# 	plot!(sol)
# end
# plot!(legend=false)


# # Controller conditions
# stops = collect(0.0:5.0:100.0)
#
#
#
#
#
# function condition(u,t,integrator)
# 	#integrator.t in tstops
# 	t in (stops)
# 	#(mod(integrator.t,durOn+durOff)-durOn)*mod(integrator.t,durOn+durOff) # condition when dynamics change
# end
# function affect!(integrator) # effect that occurs when the condition is met
# 	integrator.u[3] = 100
# end
# cb = ContinuousCallback(condition,affect!)
# prob = DiscreteProblem([init_pop; 0; 10],(0.0,time),rates)
# jump_prob = JumpProblem(prob,Direct(),LV_model)
# sol = solve(jump_prob,FunctionMap(),callback=cb,tstops=stops)
# plot(sol)
