# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
using Statistics, OrdinaryDiffEq, DiffEqJump, DiffEqBase, DifferentialEquations, Catalyst

global γ = 1e-5 # inverse carrying capacity
global κ = 1000 # competition parameter κ

global num_sims = 10 # number of realisations
global time = 100 # length of a simulation
global X₀ = 1e4 # initial number of susceptible, X

# Parameters for step lengths and periods
global step_length = 5
global num_steps = 20
global period = step_length*num_steps


# Parameters
a₁ = 4.75#0.3 # death rate of susceptible from antibiotic
a₂ = 0.1*a₁ # death rate of resistant from antibiotic
b₁ = 1 #0.2 # birth rate
d = 0.25 #0.05 # death rate
c = (1-0.88)*(b₁-d)#0.02 # cost of resistance
b₂ = b₁-c
ω = γ/(1+exp(-κ*c)) # relative fitness of susceptible in the antibiotic-free regime
ωₐ = γ/(1+exp(-κ*(c-0.9*a₁))) # relative fitness of susceptible in the antibiotic regime
Δω = ω-ωₐ
γ₁₁ = γ
γ₁₂ = γ/10
γ₂₁ = γ*10
ϝ₁₂ = γ*10
ϝ₂₁ = γ/10
γ₂₂ = γ*10
μ = 1e-4#2e-5 # mutation rate
rates = (a₁, a₂, b₁, b₂, d, γ₁₁, γ₁₂, γ₂₁, ϝ₁₂, ϝ₂₁, γ₂₂, μ)

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	b₁, X → 2X # X birth
	b₂, Y → 2Y # Y birth
    d, X → 0 # X death
	d, Y → 0 # Y death
	γ₁₁, 2X → X # self competition of X
	Ā*γ₁₂ + (1-Ā)*ϝ₁₂, X+Y → Y # competition to X from Y
	Ā*γ₂₁ + (1-Ā)*ϝ₂₁, X+Y → X # competition to Y from X
	γ₂₂, 2Y → Y # self competition of Y
	a₁*Ā, X → 0 # death by antibiotic
	a₂*Ā, Y → 0 # death by antibiotic, α₂ < α₁
	(10*Ā + (1-Ā))*μ, X → Y # mutation X to Y
	μ, Y → X # mutation Y to X
	0, Ā → 0
end a₁ a₂ b₁ b₂ d γ₁₁ γ₁₂ γ₂₁ ϝ₁₂ ϝ₂₁ γ₂₂ μ

# Solve system for various protocols averaged 100 times
output_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
time_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
time_std = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
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
			(u[1]+u[2] == 0) | ((u[2] > 50) | (u[1] > 5000))
		end
		function affectOn!(integrator) # effect that occurs when the condition is met
			integrator.u[3] = 1
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
			prob = DiscreteProblem([X₀; 0; 1],(0.0,time),rates)
			jump_prob = JumpProblem(LV_model,prob,Direct())
			sol = solve(jump_prob,FunctionMap(),callback=cbs,tstops=protocol)
			result[sim] = 1 - min(sol[1,end] + sol[2,end], 1.0)
			#result[sim] = floor(0.5*(sign(time - sol.t[end])+1)) #sum(sol[1:2,end])
			result_time[sim] = 1-sol.t[end]/time
		end
		output_mean[count,:] = [durOn durOff mean(result)]
		time_mean[count,:] = [durOn durOff mean(result_time)]
		time_std[count,:] = [durOn durOff std(result_time)]
		global count += 1
	end
end

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
gamma = γ
kappa = κ
using RCall
@rput output_mean time_mean time_std gamma kappa
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output_mean <- as.data.frame(output_mean)
time_mean <- as.data.frame(time_mean)
time_std <- as.data.frame(time_std)

p <- ggplot() + theme(strip.background=element_blank(),legend.justification = c(1, 0),
  legend.position = c(1, 0.6),plot.margin = unit(c(0, 0, 1.25, 0), "cm"),legend.title=element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top"))

p_mean <- p + geom_raster(data=output_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = "On", y = "Off") + ggtitle("Mean")

pt_mean <- p + geom_raster(data=time_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = "On", y = "Off") + ggtitle("Mean time")
pt_std <- p + geom_raster(data=time_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = "On", y = "Off") + ggtitle("SD time")

save_plot(plot_grid(p_mean,pt_mean,pt_std,ncol=3),filename=paste("/Users/brycemorsky/Desktop/",paste("test gamma",gamma,"kappa",kappa,sep="_"),".png", sep=""),base_height = 4,base_width = 12)
"""
