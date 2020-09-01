# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
using DiffEqBiological, DifferentialEquations, Statistics

# Initial conditions
control = 1.0 # amount of antibiotic when on
num_sims = 100 # number of sims to average over
time = 2000 # length of process
X₀ = 1e3 # initial number of susceptible, n₁

# Parameters
β₁ = 0.2 # susceptible birth rate
δ₁ = 0.18 # susceptible death rate
β₂ = 0.1 # resistant birth rate
δ₂ = 0.09 # resistant death rate
γ₁₁ = 1e-4 # susceptible death rate from other susceptible
γ₂₂ = 1e-4 # resistant death rate from other resistant
α₁ = 0.03 # susceptible death rate from antibiotic
α₂ = 0.0 # resistant death rate from antibiotic
μ = 2e-5 # mutation rate
η = 0.0 # 4e-3 # pasmid transfer rate

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
    β₁, n₁ → 2n₁ # n₁ birth
    δ₁, n₁ → 0 # n₁ death
	β₂, n₂ → 2n₂ # n₂ birth
	δ₂, n₂ → 0 # n₂ death
	γ₁₁, n₁ + n₁ → n₁ # self competition of n₁
	(a*κ + (1-a)/κ)*γ₁₁, n₁ + n₂ → n₂ # competition to n₁ from n₂
	(a/κ + (1-a)*κ)*γ₂₂, n₂ + n₁ → n₁ # competition to n₂ from n₁
	γ₂₂, n₂ + n₂ → n₂ # self competition of n₂
	α₁, a + n₁ → a # death by antibiotic
	α₂, a + n₂ → a # death by antibiotic, α₂ < α₁
	μ, n₁ → n₂ # mutation
	μ, n₂ → n₁ # mutation
	10*μ, a + n₁ → a + n₂ # mutation due to stress from the antibiotic
	η, n₁ + n₂ → 2n₂ # plasmid (resistance) transfer
end β₁ δ₁ β₂ δ₂ γ₁₁ κ γ₂₂ α₁ α₂ μ η

# Solve system for various protocols
function runSim()
	step_length = 5
	num_steps = 20
	period = step_length*num_steps

	output_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
	time_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)

	count = 1
	for durOn = step_length:step_length:step_length*(num_steps-1)
		for durOff = step_length:step_length:step_length*num_steps-durOn
			result = zeros(num_sims)
			result_time = zeros(num_sims)
			# Controller conditions
			protocolOn = collect(0:durOn+durOff:time)
			protocolOff = collect(durOn:durOn+durOff:time)
			# protocol = merge(protocolOn,protocolOff)
			protocol = collect(0:1:time)
			function conditionOn(u,t,integrator)
				t in (protocolOn)
			end
			function conditionOff(u,t,integrator)
				t in (protocolOff)
			end
			function conditionStop(u,t,integrator)
				(u[1]+u[2] == 0) | (u[2] > 50)
			end
			function affectOn!(integrator) # effect that occurs when the condition is met
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
			for sim = 1:num_sims
				prob = DiscreteProblem([X₀; 0; control],(0.0,time),rates)
				jump_prob = JumpProblem(prob,Direct(),LV_model)
				sol = solve(jump_prob,FunctionMap(),callback=cbs,tstops=protocol)
				result[sim] = 1 - min(sol[1,end] + sol[2,end], 1.0)
				result_time[sim] = 1-sol.t[end]/time
			end
			output_mean[count,:] = [durOn durOff mean(result)]
			time_mean[count,:] = [durOn durOff mean(result_time)]
			count += 1
		end
	end
return output_mean,time_mean
end

# Sweep the protocol
best_protocol = zeros(20,4)
best_protocol_time = zeros(20,4)

global count = 1
for κ = 1:1:20
	global rates = (β₁, δ₁, β₂, δ₂, γ₁₁, κ, γ₂₂, α₁, α₂, μ, η)
	output_mean,time_mean = runSim()
	best = findmax(output_mean[:,3])[2]
	global best_protocol[count,:] = [κ output_mean[best,:]']
	global best_protocol_time[count,:] = [κ time_mean[best,:]']
	global count +=1
end

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
using RCall
@rput best_protocol best_protocol_time
R"""

best_protocol <- as.data.frame(best_protocol)
names(best_protocol) <- c("Competition","On","Off","Success")
best_protocol_time <- as.data.frame(best_protocol_time)
names(best_protocol_time) <- c("Competition","On","Off","Time")

#save(best_protocol,file="/home/bmorsky/antibiotic/output/best_protocol.Rda")
#save(best_protocol_time,file="/home/bmorsky/antibiotic/output/best_protocol_time.Rda")

save(best_protocol,file="~/Desktop/best_protocol.Rda")
save(best_protocol_time,file="~/Desktop/best_protocol_time.Rda")
"""
