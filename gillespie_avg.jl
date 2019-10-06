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
control = 10 # initial input rate of bacteria
time = 200 # length of process
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

# Solve system for various protocols averaged 100 times
function runSim()
	output_mean = zeros(21*21,3)
	output_std = zeros(21*21,3)
	count = 1
	for durOn = 0:20
		for durOff = 0:20
			result = zeros(20)
			# Controller conditions
			function condition(t,u,integrator)
				(mod(integrator.t,durOn+durOff)-durOn)*mod(integrator.t,durOn+durOff) # condition when dynamics change
			end
			function affect!(t,u,integrator) # effect that occurs when the condition is met
				if integrator.control == 10
					integrator.control = 0
				else
					integrator.control = 10
				end
			end
			cb = ContinuousCallback(condition,affect!)
			for sim = 1:20
				prob = DiscreteProblem([init_pop; 0; 10],(0.0,time),rates)
				jump_prob = JumpProblem(prob,Direct(),LV_model)
				sol = solve(jump_prob,FunctionMap(),callback=cb,save_everystep=false,save_start=false)
				result[sim] = sum(sol[1:2,end])
				#output[count,:] += [sum(sol[1:2,end]) durOn durOff]
			end
			output_mean[count,:] = [durOn durOff mean(result)]
			output_std[count,:] = [durOn durOff std(result)]
			count += 1
		end
	end
return output_mean,output_std
end

output_mean,output_std = runSim()

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
using RCall
@rput output_mean output_std
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output_mean <- as.data.frame(output_mean)
output_std <- as.data.frame(output_std)

p <- ggplot() + theme(strip.background=element_blank()) + theme(plot.margin = unit(c(0, 0, 1.25, 0), "cm"),legend.title=element_blank()) + coord_equal() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top"))

p_mean <- p + geom_raster(data=output_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1000)) + labs(x = "On", y = "Off") + ggtitle("Mean")
p_std <- p + geom_raster(data=output_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1000)) + labs(x = "On", y = "Off") + ggtitle("SD")

save_plot(plot_grid(p_mean,p_std,ncol=2),filename="~/Documents/Notre Dame/ND paper 2/Code/heatmap.png",base_height = 4,base_width = 8)
"""
