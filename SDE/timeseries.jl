# code to generate time series (Figure 1)
using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, StochasticDiffEq

# Parameters
const ξ = 0 # stochasticity ξ
α = 1.6 # competition parameter α
const bx = ξ+1.4#0.35 #0.3# susceptible birth rate
const by = ξ+1.2#0.3 #0.25# resistant birth rate
const d = ξ+0.4#0.1 #0.05# death rate
const T = 336 #1344 # length of a simulation
const m = 10 # mutant stress
const κ = 4 # competition parameter κ
const γ = 1e-10 #1e-5
const γxx = γ # susceptible death rate from other susceptible
const γxy = γ/κ # susceptible death rate from resistant
const γyx = γ*κ # resistant death rate from susceptible
const γyy = γ # resistant death rate from other resistant
const μ = 1e-9# 1e-5 # mutation rate μ
const λ = 0.0

const X₀ = 1000000000
const A₀ = 10000
rates = [bx, by, d, γxx, γxy, γyx, γyy, α, μ, λ, m]

protocol = collect(0:0.05:T)

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	bx, X → 2X # X birth
	by, Y → 2Y # Y birth
    d, X → 0 # X death
	d, Y → 0 # Y death
	γxx*X, X → 0 # self competition of X
	γxy*Y, X → 0 # competition to X from Y
 	γyx*X, Y → 0 # competition to Y from X
	γyy*Y, Y → 0 # self competition of Y
	0.0001*A*α, X → 0 # death by antibiotic
	0.0001*A*α/10, Y → 0 # death by antibiotic
	(0.0001*(m-1)*A+1)*μ, X → Y # mutation
	μ, Y → X # mutation
	λ, A → 0
end bx by d γxx γxy γyx γyy α μ λ m

# Switching protocols
# Controller conditions
durOn = 48
durOff = 48
# Controller conditions
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
cbsPulse = CallbackSet(cbOn,cbOff,cbNonnegX,cbNonnegY,cbEnd)
cbsConst = CallbackSet(cbNonnegX,cbNonnegY,cbEnd)

prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),rates)
sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsPulse)
output1 = [sol[1,:] sol.t[:] zeros(length(sol[1,:]),1); sol[2,:] sol.t[:] ones(length(sol[2,:]),1);] #output

prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),rates)
sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsConst)
output2 = [sol[1,:] sol.t[:] zeros(length(sol[1,:]),1); sol[2,:] sol.t[:] ones(length(sol[2,:]),1);] #output
minimum(output2)
output1[output1[:,1] .< 0, 1] .= 0

#R code to generate figure
using RCall
@rput output1 output2 time
R"""
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_cowplot());

output1 <- as.data.frame(output1)
names(output1) <- c("Number","Time","Type")
output2 <- as.data.frame(output2)
names(output2) <- c("Number","Time","Type")

output1 <- ggplot() + geom_line(data=output1,aes(x=Time,y=Number,group=Type),size=1) +
ggtitle("On/off protocol") + labs(x = "Time (days)", y = bquote("Bacterial load"~(10^9))) +
scale_color_manual(values=c("blue","red"),labels = c("Wild-type", "Mutant")) +
aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,1)) +
scale_x_continuous(breaks=c(0,48,96,144,192,240,288,336),labels=c(0,2,4,6,8,10,12,14),expand = c(0.005,0)) +
scale_y_continuous(breaks=c(0,2000000000,4000000000,6000000000,8000000000,10000000000),labels=c(0,2,4,6,8,10),limits = c(0, 10001000000),expand = c(0.02,0))

output2 <- ggplot() + geom_line(data=output2,aes(x=Time,y=Number,group=Type),size=1) +
ggtitle("Constant application") + labs(x = "Time (days)", y = bquote("Bacterial load"~(10^9))) +
scale_color_manual(values=c("blue","red"),labels = c("Wild-type", "Mutant")) +
aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,1)) +
scale_x_continuous(breaks=c(0,48,96,144,192,240,288,336),labels=c(0,2,4,6,8,10,12,14),expand = c(0.005,0)) +
scale_y_continuous(breaks=c(0,2000000000,4000000000,6000000000,8000000000,10000000000),labels=c(0,2,4,6,8,10),limits = c(0, 10000000000),expand = c(0.02,0))

save_plot(output1,filename="ts24hrperiod.png",base_height = 2.5,base_width = 10)
save_plot(output2,filename="tsOn.png",base_height = 2.5,base_width = 10)
"""
