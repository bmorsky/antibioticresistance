# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
using DiffEqBiological, DifferentialEquations, Statistics


num_sims = 10

# Initial conditions
control = 1.0 # amount of antibiotic when on
time = 1000 # length of process
X₀ = 1e3 # initial number of susceptible, n₁

# Parameters
β₁ = 0.2 # susceptible birth rate
β₂ = 0.15 # resistant birth rate
δ₁ = 0.05 # susceptible death rate
δ₂ = 0.05 # resistant death rate
γ₁₁ = 1e-5 # susceptible death rate from other susceptible
γ₁₂ = 1e-5 # susceptible death rate from resistant
γ₂₁ = 1e-5 # resistant death rate from susceptible
γ₂₂ = 1e-5 # resistant death rate from other resistant
α₁ = 0.3 # susceptible death rate from antibiotic
α₂ = 0.03 # resistant death rate from antibiotic
μ = 2e-5 # mutation rate
η = 0.0 # 4e-3 # pasmid transfer rate
rates = (β₁, β₂, δ₁, δ₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, η)

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	β₁, n₁ → 2n₁ # n₁ birth
	β₂, n₂ → 2n₂ # n₂ birth
    δ₁, n₁ → 0 # n₁ death
	δ₂, n₂ → 0 # n₂ death
	γ₁₁, n₁ + n₁ → n₁ # self competition of n₁
	γ₁₂, n₁ + n₂ → n₂ # competition to n₁ from n₂
	γ₂₁, n₂ + n₁ → n₁ # competition to n₂ from n₁
	γ₂₂, n₂ + n₂ → n₂ # self competition of n₂
	control*α₁, n₁ → 0 # death by antibiotic
	control*α₂, n₂ → 0 # death by antibiotic, α₂ < α₁
	μ, n₁ → n₂ # mutation
	(control*10+(1-control))*μ, n₂ → n₁ # mutation
	#a*10*μ, n₁ → n₂ # mutation due to stress from the antibiotic
	η, n₁ + n₂ → 2n₂ # plasmid (resistance) transfer
	0, control → 0
end β₁ β₂ δ₁ δ₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ η

# Solve system for various protocols averaged 100 times
step_length = 10
num_steps = 10
period = step_length*num_steps

output_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
output_std = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
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
			#result[sim] = floor(0.5*(sign(time - sol.t[end])+1)) #sum(sol[1:2,end])
			result_time[sim] = 1-sol.t[end]/time
		end
		output_mean[count,:] = [durOn durOff mean(result)]
		output_std[count,:] = [durOn durOff std(result)]
		time_mean[count,:] = [durOn durOff mean(result_time)]
		time_std[count,:] = [durOn durOff std(result_time)]
		global count += 1
	end
end

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
using RCall
@rput output_mean output_std time_mean time_std
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output_mean <- as.data.frame(output_mean)
output_std <- as.data.frame(output_std)
time_mean <- as.data.frame(time_mean)
time_std <- as.data.frame(time_std)

p <- ggplot() + theme(strip.background=element_blank(),legend.justification = c(1, 0),
  legend.position = c(1, 0.6),plot.margin = unit(c(0, 0, 1.25, 0), "cm"),legend.title=element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top"))

p_mean <- p + geom_raster(data=output_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = "On", y = "Off") + ggtitle("Mean")
p_std <- p + geom_raster(data=output_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = "On", y = "Off") + ggtitle("SD")

pt_mean <- p + geom_raster(data=time_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = "On", y = "Off") + ggtitle("Mean time")
pt_std <- p + geom_raster(data=time_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = "On", y = "Off") + ggtitle("SD time")

save_plot(plot_grid(p_mean,p_std,pt_mean,pt_std,ncol=2),filename="~/Desktop/heatmap2.png",base_height = 8,base_width = 8)
#save_plot(plot_grid(p_mean,p_std,pt_mean,pt_std,ncol=2),filename="/home/bmorsky/antibiotic/heatmap2.png",base_height = 8,base_width = 8)
"""
