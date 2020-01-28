# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol.
using DiffEqBiological, DifferentialEquations, Statistics

num_sims = 10
control = 1 # initial input rate of bacteria
time = 1000 # length of process
# X₀ = 1e4 # initial number of susceptible, n₁

# Stochastic LV chemical reaction system
# Parameters, times are per hour

# Parameters
β₁ = 0.02 # susceptible birth rate
δ₁ = 0.014 # susceptible death rate
β₂ = 0.01 # resistant birth rate
δ₂ = 0.005 # resistant death rate
γ₁₁ = 1e-6 # susceptible death rate from other susceptible
γ₁₂ = 0.88e-6 # susceptible death rate from resistant
γ₂₁ = 1.36e-6 # resistant death rate from susceptible
γ₂₂ = 1e-6 # resistant death rate from other resistant
α₁ = 0.008 # susceptible death rate from antibiotic
α₂ = 0.0015 # resistant death rate from antibiotic
μ = 0.10e-6 # mutation rate
μₐ = 10*μ # mutation rate cause by antibiotic
η = 0.0 #4e-3 # pasmid transfer rate
X₀ = 5e3 #initial population size
final_time = 5000

# b₁ = 0.02 # susceptible birth rate
# d₁ = 0.014 # susceptible death rate
# b₂ = 0.01 # resistant birth rate
# d₂ = 0.005 # resistant death rate
# γ₁₁ = 1e-6 # susceptible death rate from other susceptible
# γ₁₂ = 1e-6 # susceptible death rate from resistant
# γ₂₁ = 1e-6 # resistant death rate from susceptible
# γ₂₂ = 1e-6 # resistant death rate from other resistant
# α₁ = 0.0008 # susceptible death rate from antibiotic
# α₂ = 0.00015 # resistant death rate from antibiotic
# μ = 1e-6 # mutation rate
# μₐ = 10*μ # mutation rate cause by antibiotic
# η = 0 #4e-3 # pasmid transfer rate
# init_pop = 1e4 #initial population size

# Initial conditions
antibiotic_amount = 1.0 # amount of antibiotic when on
# durOn = 700 # durtation the antibiotic is on
# durOff = 1500 # duration the antibiotic is off

# b₁ = 0.04
# d₁ = 0.03
# b₂ = 0.03
# d₂ = 0.02
# γ₁₁ = 1e-4
# γ₁₂ = 0.88e-4
# γ₂₁ = 1.136e-4
# γ₂₂ = 1e-4
# α₁ = 0.0012
# α₂ = 0
# μ = 1.5e-8
# ζ = 0.0
# η = 0 #4e-3
# controlrate = 0
rates = (β₁, δ₁, β₂, δ₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, μₐ, η)
LV_model = @reaction_network begin
    β₁, n₁ → 2n₁ # n₁ birth
    δ₁, n₁ → 0 # n₁ death
	β₂, n₂ → 2n₂ # n₂ birth
	δ₂, n₂ → 0 # n₂ death
	γ₁₁, n₁ + n₁ → 0 # self competition of n₁
	γ₁₂, n₁ + n₂ → 0 # competition to n₁ from n₂
	γ₂₁, n₂ + n₁ → 0 # competition to n₂ from n₁
	γ₂₂, n₂ + n₂ → 0 # self competition of n₂
	α₁, a + n₁ → a # death by antibiotic
	α₂, a + n₂ → a # death by antibiotic, α₂ < α₁
	μ, n₁ → n₂ # mutation
	μ, n₂ → n₁ # mutation
	μₐ, a + n₁ → a + n₂ # mutation due to stress from the antibiotic
	η, n₁ + n₂ → 2n₂ # plasmid (resistance) transfer
end β₁ δ₁ β₂ δ₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ η

# Solve system for various protocols averaged 100 times
step_length = 1
num_steps = 1000
period = step_length*num_steps
function runSim()
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
			count += 1
		end
	end
return output_mean,output_std,time_mean,time_std
end
# Sweep the protocols
output_mean,output_std,time_mean,time_std = runSim()
# Sweep X₀ and μ
output_muX0_mean = zeros(10^2,3)
output_muX0_std = zeros(10^2,3)
time_muX0_mean = zeros(10^2,3)
time_muX0_std = zeros(10^2,3)
# global count = 1
# for X₀ = 500:100:1500
# 	for μ = 1e-6:1e-6:1e-5
# 		rates = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, ζ, η, controlrate)
# 		output_mean,output_std,time_mean,time_std = runSim()
# 		output_muX0_mean[count,:] = [X₀ μ minimum(output_mean[:,3])]
# 		output_muX0_std[count,:] = [X₀ μ output_std[argmin(output_mean[:,3]),3]]
# 		time_muX0_mean[count,:] = [X₀ μ minimum(time_mean[:,3])]
# 		time_muX0_std[count,:] = [X₀ μ time_std[argmin(time_mean[:,3]),3]]
# 		global count +=1
# 	end
# end

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
using RCall
@rput output_mean output_std time_mean time_std output_muX0_mean output_muX0_std time_muX0_mean time_muX0_std
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output_mean <- as.data.frame(output_mean)
output_std <- as.data.frame(output_std)
time_mean <- as.data.frame(time_mean)
time_std <- as.data.frame(time_std)
output_muX0_mean <- as.data.frame(output_muX0_mean)
output_muX0_std <- as.data.frame(output_muX0_std)
time_muX0_mean <- as.data.frame(time_muX0_mean)
time_muX0_std <- as.data.frame(time_muX0_std)

p <- ggplot() + theme(strip.background=element_blank(),legend.justification = c(1, 0),
  legend.position = c(1, 0.6),plot.margin = unit(c(0, 0, 1.25, 0), "cm"),legend.title=element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top"))

p_mean <- p + geom_raster(data=output_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = "On", y = "Off") + ggtitle("Mean")
p_std <- p + geom_raster(data=output_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = "On", y = "Off") + ggtitle("SD")

pt_mean <- p + geom_raster(data=time_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = "On", y = "Off") + ggtitle("Mean time")
pt_std <- p + geom_raster(data=time_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = "On", y = "Off") + ggtitle("SD time")


q <- ggplot() + theme(plot.margin = unit(c(0, 0, 1.25, 0), "cm"),legend.title=element_blank()) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top"))

q_muX0_mean <- q + geom_raster(data=output_muX0_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = expression(X[0]), y = expression(mu)) + ggtitle("Mean")
q_muX0_std <- q + geom_raster(data=output_muX0_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = expression(X[0]), y = expression(mu)) + ggtitle("SD")
q_time_muX0_mean <- q + geom_raster(data=time_muX0_mean,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis",limits = c(0,1)) + labs(x = expression(X[0]), y = expression(mu)) + ggtitle("Mean time")
q_time_muX0_std <- q + geom_raster(data=time_muX0_std,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="viridis") + labs(x = expression(X[0]), y = expression(mu)) + ggtitle("SD time")

save_plot(plot_grid(p_mean,p_std,pt_mean,pt_std,ncol=2),filename="/home/bmorsky/antibiotic/heatmap.png",base_height = 8,base_width = 8)
save_plot(plot_grid(q_muX0_mean,q_muX0_std,q_time_muX0_mean,q_time_muX0_std,ncol=2),filename="/home/bmorsky/antibiotic/heatmap_muX0.png",base_height = 8,base_width = 8)
"""









# Solve system for various protocols averaged 100 times
# function runSim()
# 	output_mean = zeros(21*21,3)
# 	output_std = zeros(21*21,3)
# 	count = 1
# 	for durOn = 0:20
# 		for durOff = 0:20
# 			result = zeros(20)
# 			# Controller conditions
# 			function condition(t,u,integrator)
# 				(mod(integrator.t,durOn+durOff)-durOn)*mod(integrator.t,durOn+durOff) # condition when dynamics change
# 			end
# 			function affect!(t,u,integrator) # effect that occurs when the condition is met
# 				if integrator.control == 10
# 					integrator.control = 0
# 				else
# 					integrator.control = 10
# 				end
# 			end
# 			cb = ContinuousCallback(condition,affect!)
# 			for sim = 1:20
# 				prob = DiscreteProblem([init_pop; 0; 10],(0.0,time),rates)
# 				jump_prob = JumpProblem(prob,Direct(),LV_model)
# 				sol = solve(jump_prob,FunctionMap(),callback=cb,save_everystep=false,save_start=false)
# 				result[sim] = sum(sol[1:2,end])
# 				#output[count,:] += [sum(sol[1:2,end]) durOn durOff]
# 			end
# 			output_mean[count,:] = [durOn durOff mean(result)]
# 			output_std[count,:] = [durOn durOff std(result)]
# 			count += 1
# 		end
# 	end
# return output_mean,output_std
# end
