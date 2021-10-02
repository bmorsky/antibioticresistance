using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, OrdinaryDiffEq, Statistics, StochasticDiffEq

# Parameters
const α = 1.6 # death rate from antibiotic
const b = 1.4 # birth rate
const d = 0.4 # death rate
const m = 10 # mutant stress
const num_sims = 200 # number of realisations
const T = 336 # length of a time in hours (2 weeks)
const κ = 4 # competition parameter κ
const λ = 0.0 # antibiotic degredation rate
const μ = 1e-9 # mutation rate μ

# Initial conditions
const X₀ = 1000000000
const A₀ = 10000

# Output
output = zeros(2*num_sims*100,3)
# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	bx, X → 2X # X birth
	by, Y → 2Y # Y birth
    d, X → 0 # X death
	d, Y → 0 # Y death
	γ*X, X → 0 # self competition of X
	(γ/κ)*Y, X → 0 # competition to X from Y
 	κ*γ*X, Y → 0 # competition to Y from X
	γ*Y, Y → 0 # self competition of Y
	0.0001*A*α, X → 0 # death by antibiotic
	0.00001*A*α, Y → 0 # death by antibiotic
	(0.0001*(m-1)*A+1)*μ, X → Y # mutation
	μ, Y → X # mutation
	λ, A → 0
end bx by d γ κ α μ λ m

# Switching protocols
protocol = collect(0:0.05:T) # time stops for callbacks
const durOn = 2 # duration on
const durOff = 2 # duration off
const protocolOn = collect(0:durOn+durOff:T)
const protocolOff = collect(durOn:durOn+durOff:T)
conditionOn(u,t,integrator) = t in (protocolOn)
conditionOff(u,t,integrator) = t in (protocolOff)
conditionNonnegX(u,t,integrator) = u[1]<1
conditionNonnegY(u,t,integrator) = u[2]<1
conditionEnd(u,t,integrator) = (u[1]<1 && u[2]<1)
function affectOn!(integrator) # effect that occurs when the condition is met
	#integrator.p[10] = 0
	integrator.u[3] = A₀
	reset_aggregated_jumps!(integrator)
end
function affectOff!(integrator) # effect that occurs when the condition is met
	integrator.u[3] = 0
	reset_aggregated_jumps!(integrator)
end
function affectNonnegX!(integrator) # effect that occurs when the condition is met
	integrator.u[1] = 0
	reset_aggregated_jumps!(integrator)
end
function affectNonnegY!(integrator) # effect that occurs when the condition is met
	integrator.u[2] = 0
	reset_aggregated_jumps!(integrator)
end
affectEnd!(integrator) = terminate!(integrator) # effect that occurs when the condition is met
cbOn = DiscreteCallback(conditionOn,affectOn!,save_positions = (false,false))
cbOff = DiscreteCallback(conditionOff,affectOff!,save_positions = (false,false))
cbNonnegX = DiscreteCallback(conditionNonnegX,affectNonnegX!,save_positions = (false,false))
cbNonnegY = DiscreteCallback(conditionNonnegY,affectNonnegY!,save_positions = (false,false))
cbEnd = DiscreteCallback(conditionEnd,affectEnd!,save_positions = (false,false))
cbsPulse = CallbackSet(cbOn,cbOff,cbNonnegX,cbNonnegY,cbEnd,save_positions = (false,false))
cbsConst = CallbackSet(cbNonnegX,cbNonnegY,cbEnd,save_positions = (false,false))

function runPulse(γ)
	prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),[bx, by, d, γ, κ, α, μ, λ, m])
	sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsPulse)
	return sum(sol[1:2,:]/(10*T+1))
end

function runConst(γ)
	prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),[bx, by, d, γ, κ, α, μ, λ, m])
	sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsConst)
	return sum(sol[1:2,:]/(10*T+1))
end

function ensemPulse!(γ)
	for sim = 1:num_sims
		global output[count,:] = [γ runPulse(γ) 1];
		global count = count + 1;
	end
end

function ensemConst!(γ)
	for sim = 1:num_sims
 		global output[count,:] = [γ runConst(γ) 0];
		global count = count + 1;
	end
end

function runSim!()
	for γ = 1e-11:1e-11:1e-9
		ensemPulse!(γ)
		ensemConst!(γ)
	end
end

# Run simulations
global count = 1
runSim!()

#R code to generate gamma vs relative average bacterial load
using RCall
@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridisLite)
library(viridis)

output <- as.data.frame(output)
names(output) <- c("Parameter","Result","Type")
save(output,file="gamma.Rda")

output <- ggplot(data=output,aes(x=Parameter,y=Result,group=Type)) + geom_point(size=1) +
stat_summary(fun=mean, geom="line") +
labs(x = bquote(gamma~(10^{-11})), y = bquote("Bacterial load"~(10^{10}))) +
scale_color_manual(values=c("red","blue"),labels = c("Constant", "Pulsed")) +
aes(color = factor(Type)) +
theme(legend.title=element_blank(),legend.position="right",legend.direction='vertical') +
scale_x_continuous(breaks=c(0.00000000001,0.0000000002,0.0000000004,0.0000000006,0.0000000008,0.000000001),labels=c(1,20,40,60,80,100),limits = c(0.00000000001,0.000000001),expand = c(0.005,0)) +
scale_y_continuous(breaks=c(0,20000000000,40000000000,60000000000),labels=c(0,2,4,6),limits = c(0,60000000000),expand = c(0.02,0))

save_plot(output,filename="gamma.png",base_height = 2.5,base_width = 8)
"""
