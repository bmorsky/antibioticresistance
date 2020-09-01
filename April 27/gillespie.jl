# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
# Packages
using DifferentialEquations, DiffEqBiological, Plots

# Initial conditions
control = 1.0 # amount of antibiotic when on
num_sims = 100 # number of sims to average over
time = 2000 # length of process
X₀ = 1e3 # initial number of susceptible, n₁

# Parameters
α₁ = 0.03 # susceptible death rate from antibiotic
α₂ = 0.0 # resistant death rate from antibiotic
β₁ = 0.2 # susceptible birth rate
β₂ = 0.1 # resistant birth rate
γ₁₁ = 1e-5 # susceptible death rate from other susceptible
γ₂₂ = 1e-5 # resistant death rate from other resistant
δ₁ = 0.18 # susceptible death rate
δ₂ = 0.09 # resistant death rate
κ = 4#2
μ = 2e-5 # mutation rate
η = 0.0 # 4e-3 # pasmid transfer rate

durOn = 50 # durtation the antibiotic is on
durOff = 25 # duration the antibiotic is off

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
rates = (β₁, δ₁, β₂, δ₂, γ₁₁, κ, γ₂₂, α₁, α₂, μ, η)

# Controller conditions
protocolOn = collect(0:durOn+durOff:time)
protocolOff = collect(durOn:durOn+durOff:time)
protocol = merge(protocolOn,protocolOff)
protocol = collect(0:1:time)
# Conditions to change the integrator
function conditionOn(u,t,integrator) # turn the antibiotic on
	t in (protocolOn)
end
function conditionOff(u,t,integrator) # turn the antibiotic off
	#false
	t in (protocolOff)
end
function conditionStop(u,t,integrator) # when the bacteria are eliminated
	u[1]+u[2] == 0
end
# Effects that occurs when the condition is met
function affectOn!(integrator)
	integrator.u[3] = control # turn on the antibiotic
end
function affectOff!(integrator)
	integrator.u[3] = 0 # turn off the antibiotic
end
function affectStop!(integrator)
	terminate!(integrator) # end solver
end
cbOn = DiscreteCallback(conditionOn,affectOn!)
cbOff = DiscreteCallback(conditionOff,affectOff!)
cbStop = DiscreteCallback(conditionStop,affectStop!)
cbs = CallbackSet(cbOn,cbOff,cbStop)
prob = DiscreteProblem([X₀; 0; control],(0.0,time),rates)
jump_prob = JumpProblem(prob,Direct(),LV_model)
sol1 = solve(jump_prob,FunctionMap(),callback=cbs,tstops=protocol)
sol2 = solve(jump_prob,FunctionMap(),callback=cbStop,tstops=protocol)

# Plot the solution. The blue line is the susceptibles, the red resistant
# mutants, and the green the antibiotic.
#plot(sol)
#for plotting in R
output1 = [sol1[1,:] sol1.t[:] zeros(length(sol1[1,:]),1); sol1[2,:] sol1.t[:] ones(length(sol1[2,:]),1)] #output
output2 = [sol2[1,:] sol2.t[:] zeros(length(sol2[1,:]),1); sol2[2,:] sol2.t[:] ones(length(sol2[2,:]),1)] #output

#R code to generate figure
using RCall
@rput output1 output2;
R"""
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot());

output1 <- as.data.frame(output1)
names(output1) <- c("Number","Time","Type")
output2 <- as.data.frame(output2)
names(output2) <- c("Number","Time","Type")

ts <- ggplot()

output1 <- ts + geom_line(data=output1,aes(x=Time,y=Number,group=Type),size=1) + ggtitle("On/off protocol")+ scale_color_manual(values=c("blue","red"),labels = c("Susceptible", "Resistant")) + aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,.55)) + labs(x = "Time", y = "Bacterial load") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(limits = c(0,2050), expand = c(0,0))
output2 <- ts + geom_line(data=output2,aes(x=Time,y=Number,group=Type),size=1) + ggtitle("Constant application")+ scale_color_manual(values=c("blue","red"),labels = c("Susceptible", "Resistant")) + aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,.55)) + labs(x = "Time", y = "Bacterial load") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(limits = c(0,2050), expand = c(0,0))

plot_out <- plot_grid(output1,output2,labels=letters[1:2],ncol=2)
save(best_protocol,file="~/Desktop/gillespie_ex.Rda")
save_plot(plot_out,filename="~/Documents/Notre Dame/ND paper 2/Code/ts_protocol_vs_constant_stochastic.png",base_height = 5,base_width = 10)
"""


# # Run multiple realizations. Solve system num_sims times and plot the results.
# num_sims = 30
# for sim = 1:num_sims
# 	prob = DiscreteProblem([init_pop; 0; 10],(0.0,time),rates)
# 	jump_prob = JumpProblem(prob,Direct(),LV_model)
# 	sol = solve(jump_prob,FunctionMap(),callback=cb)
# 	plot!(sol)
# end
# plot!(legend=false)


# Parameters
# b₁ = 0.005 # susceptible birth rate
# d₁ = 0.001 # susceptible death rate
# b₂ = 0.005 # resistant birth rate
# d₂ = 0.004 # resistant death rate
# γ₁₁ = 1e-5 # susceptible death rate from other susceptible
# γ₁₂ = 1e-5 # susceptible death rate from resistant
# γ₂₁ = 1e-5 # resistant death rate from susceptible
# γ₂₂ = 1e-5 # resistant death rate from other resistant
# α₁ = 0.0001 # susceptible death rate from antibiotic
# α₂ = 0 # resistant death rate from antibiotic
# μ = 0.10e-6 # mutation rate
# μₐ = 3*μ # mutation rate cause by antibiotic
# η = 0 #4e-3 # pasmid transfer rate

# b₁ = 0.2 # susceptible birth rate
# d₁ = 0.18 # susceptible death rate
# b₂ = 0.192 # resistant birth rate
# d₂ = 0.18 # resistant death rate
# γ₁₁ = 1e-5 # susceptible death rate from other susceptible
# γ₁₂ = 0.5e-5 # susceptible death rate from resistant
# γ₂₁ = 2e-5 # resistant death rate from susceptible
# γ₂₂ = 1e-5 # resistant death rate from other resistant
# α₁ = 0.04 # susceptible death rate from antibiotic
# α₂ = 0.0 # resistant death rate from antibiotic
# μ = 1e-5 # mutation rate
# μₐ = 10*μ # mutation rate cause by antibiotic
# η = 0 #4e-3 # pasmid transfer rate
# init_pop = 1e3 #initial population size

# Install required packages. Comment out when installed.
# using Pkg
# Pkg.add("DifferentialEquations")
# Pkg.add("DiffEqBiological")
# Pkg.add("Plots")

# Stochastic LV chemical reaction system. Susceptible bacteria n₁, resistant
# bacteria n₂, and antibiotic a. rate, reaction → product.
# LV_model = @reaction_network begin
#     b₁, n₁ → 2n₁ # n₁ birth
#     d₁, n₁ → 0 # n₁ death
# 	b₂, n₂ → 2n₂ # n₂ birth
# 	d₂, n₂ → 0 # n₂ death
# 	γ₁₁, n₁ + n₁ → 0 # self competition of n₁
# 	a*γ₂₁+(1-a)*γ₁₂, n₁ + n₂ → 0 # competition to n₁ from n₂
# 	a*γ₁₂+(1-a)*γ₂₁, n₂ + n₁ → 0 # competition to n₂ from n₁
# 	γ₂₂, n₂ + n₂ → 0 # self competition of n₂
# 	α₁, a + n₁ → a # death by antibiotic
# 	α₂, a + n₂ → a # death by antibiotic, α₂ < α₁
# 	μ, n₁ → n₂ # mutation
# 	μ, n₂ → n₁ # mutation
# 	μₐ, a + n₁ → a + n₂ # mutation due to stress from the antibiotic
# 	η, n₁ + n₂ → 2n₂ # plasmid (resistance) transfer
# end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ η

# Parameters
# b₁ = 0.04 # susceptible birth rate
# d₁ = 0.025 # susceptible death rate
# b₂ = 0.03 # resistant birth rate
# d₂ = 0.015 # resistant death rate
# γ₁₁ = 1e-5 # susceptible death rate from other susceptible
# γ₁₂ = 0.88e-5 # susceptible death rate from resistant
# γ₂₁ = 1.136e-5 # resistant death rate from susceptible
# γ₂₂ = 1e-5 # resistant death rate from other resistant
# α₁ = 0.012 # susceptible death rate from antibiotic
# α₂ = 0 # resistant death rate from antibiotic
# μ = 1e-8 # mutation rate
# μₐ = 1000*μ # mutation rate cause by antibiotic
# η = 0 #4e-3 # pasmid transfer rate
#
# # Initial conditions
# antibiotic_amount = 100.0 # amount of antibiotic when on
# durOn = 2 # durtation the antibiotic is on
# durOff = 8 # duration the antibiotic is off
# time = 1000 # length of process
# rates = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, μₐ, η)
