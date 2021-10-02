using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, OrdinaryDiffEq, Statistics, StochasticDiffEq

# Parameters
const halflife = 60
const λ = log(2)/(halflife/60)
const λvar = 0

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
const μ = 1e-9 # mutation rate μ

# Initial conditions
const X₀ = 1000000000
const A₀ = 10000

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
	0.00001*A*α, Y → 0 # death by antibiotic
	(0.0001*(m-1)*A+1)*μ, X → Y # mutation
	μ, Y → X # mutation
	λvar, A → 0
end bx by d γxx γxy γyx γyy α μ λvar λ m

# Callbacks
protocol = collect(0:0.05:T) # time stops for callbacks
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
cbNonnegX = DiscreteCallback(conditionNonnegX,affectNonnegX!,save_positions = (false,false))
cbNonnegY = DiscreteCallback(conditionNonnegY,affectNonnegY!,save_positions = (false,false))
cbEnd = DiscreteCallback(conditionEnd,affectEnd!,save_positions = (false,false))
cbsConst = CallbackSet(cbNonnegX,cbNonnegY,cbEnd,save_positions = (false,false))

# Callbacks for PAE
function affectOnPAE!(integrator) # effect that occurs when the condition is met
	integrator.p[10] = 0
	integrator.u[3] = A₀
	reset_aggregated_jumps!(integrator)
end
function affectOffPAE!(integrator) # effect that occurs when the condition is met
	integrator.p[10] = integrator.p[11]
	reset_aggregated_jumps!(integrator)
end
conditionNonnegA(u,t,integrator) = u[3]<1
function affectNonnegA!(integrator) # effect that occurs when the condition is met
	integrator.u[3] = 0
	reset_aggregated_jumps!(integrator)
end
cbNonnegA = DiscreteCallback(conditionNonnegA,affectNonnegA!,save_positions = (false,false))

# Always on
output_alwayson = 0
protocol = collect(0:1:T)
for sim = 1:num_sims
	prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),[bx, by, d, γxx, γxy, γyx, γyy, α, μ, λvar, λ, m])
	sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsConst)
	global output_alwayson += sum(sol[1:2,:]/((10*T+1)*num_sims))
end

# Switching protocols
const step_length = 1
const num_steps = 9
const period = step_length*num_steps

output_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
count = 1
for durOn = step_length:step_length:step_length*(num_steps-1)
	for durOff = step_length:step_length:step_length*num_steps-durOn
		result = 0
		# Controller conditions
		protocolOn = collect(0:durOn+durOff:T)
		protocolOff = collect(durOn:durOn+durOff:T)
		conditionOn(u,t,integrator) = t in (protocolOn)
		conditionOff(u,t,integrator) = t in (protocolOff)
		cbOn = DiscreteCallback(conditionOn,affectOn!,save_positions = (false,false))
		cbOff = DiscreteCallback(conditionOff,affectOff!,save_positions = (false,false))
		cbsPulse = CallbackSet(cbOn,cbOff,cbNonnegX,cbNonnegY,cbEnd,save_positions = (false,false))
		cbOnPAE = DiscreteCallback(conditionOn,affectOnPAE!,save_positions = (false,false))
		cbOffPAE = DiscreteCallback(conditionOff,affectOffPAE!,save_positions = (false,false))
		cbsPulsePAE = CallbackSet(cbOnPAE,cbOffPAE,cbNonnegX,cbNonnegY,cbEnd,cbNonnegA,save_positions = (false,false))

		for sim = 1:num_sims
			prob = SDEProblem(LV_model,[X₀; 0; A₀],(0.0,T),[bx, by, d, γxx, γxy, γyx, γyy, α, μ, λvar, λ, m])
			sol = solve(prob,EM(),dt=0.05,saveat=0.1,tstops=protocol,callback=cbsPulsePAE)
			result += sum(sol[1:2,:]/((10*T+1)*num_sims))
		end
		output_mean[count,:] = [durOn durOff result]
		global count += 1
	end
end

#R code to generate heatmaps for the time average of the bacterial load
output_mean[:,3] = output_mean[:,3]/output_alwayson
maxSwitching = maximum(output_mean[:,3])
output_alwayson = 1
using RCall
@rput output_mean output_alwayson maxSwitching halflife
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)
library(RColorBrewer)
library("scales")

pal <- rev(brewer.pal(11,"RdYlBu"))

output_mean <- as.data.frame(output_mean)

p_mean <- ggplot() + theme(strip.background=element_blank(),legend.justification = c(1, 0), legend.position = c(1, 0.6),plot.margin = unit(c(0, 0, 1.25, 0), "cm")) +
guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top", title="Bacterial load")) +
geom_raster(data=output_mean,aes(x=V1,y=V2,fill=V3)) +
scale_fill_gradientn(colors=pal, values = rescale(c(0,output_alwayson,maxSwitching)), limits=c(0,maxSwitching)) +
labs(x = "Time ''on'' (hours)", y = "Time ''off'' (hours)") +
scale_y_continuous(breaks=c(2,4,6,8),labels=c(2,4,6,8),expand = c(0,0)) +
scale_x_continuous(breaks=c(2,4,6,8),labels=c(2,4,6,8),expand = c(0,0)) +
ggtitle(bquote("half-life"==.(halflife)~" minutes"))

save_plot(p_mean,filename=paste(paste("halflife",halflife,sep="_"),".png", sep=""),base_height = 4,base_width = 4)
save(output_mean,file=paste(paste("halflife",halflife,sep="_"),".Rda", sep=""))
"""