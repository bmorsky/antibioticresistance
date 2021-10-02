using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, OrdinaryDiffEq, Statistics, StochasticDiffEq

# Parameters
const α = 1.6 # death rate from antibiotic
const bx = 1.4 # susceptible birth rate
const by = 1.2 # resistant birth rate
const d = 0.4 # death rate
const m = 10 # mutant stress
const num_sims = 200 # number of realisations
const T = 336 # length of a time in hours (2 weeks)
const γ = 1e-10 # interaction rate
const γxx = γ # susceptible death rate from other susceptible
const γxy = γ/κ # susceptible death rate from resistant
const γyx = γ*κ # resistant death rate from susceptible
const γyy = γ # resistant death rate from other resistant
const κ = 4 # competition parameter κ
const λ = 0.0 # antibiotic degredation rate
const μ = 1e-9 # mutation rate μ

# Initial conditions
const X₀ = 1000000000
const A₀ = 10000

# Output
output = zeros(2*num_sims*52,3)

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	bx+ξ, X → 2X # X birth
	by+ξ, Y → 2Y # Y birth
    d+ξ, X → 0 # X death
	d+ξ, Y → 0 # Y death
	γxx*X, X → 0 # self competition of X
	γxy*Y, X → 0 # competition to X from Y
 	γyx*X, Y → 0 # competition to Y from X
	γyy*Y, Y → 0 # self competition of Y
	0.0001*A*α, X → 0 # death by antibiotic
	0.00001*A*α, Y → 0 # death by antibiotic
	(0.0001*(m-1)*A+1)*μ, X → Y # mutation
	μ, Y → X # mutation
	λ, A → 0
end bx by d γxx γxy γyx γyy α μ λ m ξ

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

function runPulse(ξ)
	prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),[bx, by, d, γxx, γxy, γyx, γyy, α, μ, λ, m, ξ])
	sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsPulse)
	return sum(sol[1:2,:]/(10*T+1))
end

function runConst(ξ)
	prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),[bx, by, d, γxx, γxy, γyx, γyy, α, μ, λ, m, ξ])
	sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsConst)
	return sum(sol[1:2,:]/(10*T+1))
end

function ensemPulse!(ξ)
	for sim = 1:num_sims
		global output[count,:] = [ξ runPulse(ξ) 1];
		global count = count + 1;
	end
end

function ensemConst!(ξ)
	for sim = 1:num_sims
 		global output[count,:] = [ξ runConst(ξ) 0];
		global count = count + 1;
	end
end

function runSim!()
	for ξ = -0.2:0.2:10
		ensemPulse!(ξ)
		ensemConst!(ξ)
	end
end

# Run simulations
global count = 1
runSim!()

#R code to generate xi vs relative average bacterial load (for SI Appendix)
using RCall
@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridisLite)
library(viridis)

output <- as.data.frame(output)
names(output) <- c("Parameter","Result","Type")
save(output,file="xi.Rda")

output <- ggplot(data=output,aes(x=Parameter,y=Result,group=Type)) + geom_point(size=1) +
stat_summary(fun=mean, geom="line") +
labs(x = bquote(xi), y = bquote("Bacterial load"~(10^9))) +
scale_color_manual(values=c("red","blue"),labels = c("Constant", "Pulsed")) +
aes(color = factor(Type)) +
theme(legend.title=element_blank(),legend.position="right",legend.direction='vertical') +
scale_x_continuous(breaks=c(-0.2,2,4,6,8,10),expand = c(0.005,0)) +
scale_y_continuous(breaks=c(0,2000000000,4000000000,6000000000,8000000000),labels=c(0,2,4,6,8),limits = c(0, 8000000000),expand = c(0.02,0))

save_plot(output,filename="xi.png",base_height = 2.5,base_width = 8)
"""
