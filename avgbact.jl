using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, OrdinaryDiffEq, Statistics, ModelingToolkit

# generations where 1 generation = 1 hr

# Parameters
b₁ = ξ+0.35 #0.3# susceptible birth rate
b₂ = ξ+0.3 #0.25# resistant birth rate
d = ξ+0.1 #0.05# death rate

X₀ = 1/γ
γ₁₁ = 2*γ # susceptible death rate from other susceptible
γ₁₂ = γ/κ # susceptible death rate from resistant
γ₂₁ = γ*κ # resistant death rate from susceptible
γ₂₂ = 2*γ # resistant death rate from other resistant
rates = (b₁, b₂, d, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α, μ, m)

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	b₁, X → 2X # X birth
	b₂, Y → 2Y # Y birth
    d, X → 0 # X death
	d, Y → 0 # Y death
	γ₁₁, 2X → X # self competition of n₁
	γ₁₂, X+Y → Y # competition to n₁ from n₂
 	γ₂₁, X+Y → X # competition to n₂ from n₁
	γ₂₂, 2Y → Y # self competition of n₂
	A*α, X → 0 # death by antibiotic
	A*α/10, Y → 0 # death by antibiotic, α₂ < α₁
	μ*(1+(m-1)*A), X → Y # mutation
	μ, Y → X # mutation
	0, A → 0
end b₁ b₂ d γ₁₁ γ₁₂ γ₂₁ γ₂₂ α μ m

# Solve system for various protocols averaged 100 times
protocol = collect(0:1:fin_time)

# Always on
output_alwayson = 0
for sim = 1:num_sims
	prob = DiscreteProblem([X₀; 0; 1],(0.0,fin_time),rates)
	jump_prob = JumpProblem(LV_model,prob,Direct(),save_positions=(false,false))
	sol = solve(jump_prob,SSAStepper(),saveat=0.1)
	global output_alwayson += sum(sol[1:2,:])/((10*fin_time+1)*num_sims)
end

# Switching protocols
output_mean = zeros(convert(Int64,num_steps*(num_steps-1)/2),3)
count = 1
for durOn = step_length:step_length:step_length*(num_steps-1)
	for durOff = step_length:step_length:step_length*num_steps-durOn
		result = 0
		# Controller conditions
		protocolOn = collect(0:durOn+durOff:fin_time)
		protocolOff = collect(durOn:durOn+durOff:fin_time)
		# protocol = merge(protocolOn,protocolOff)
		conditionOn(u,t,integrator) = t in (protocolOn)
		conditionOff(u,t,integrator) = t in (protocolOff)

		function affectOn!(integrator) # effect that occurs when the condition is met
			integrator.u[3] = 1
		end
		function affectOff!(integrator) # effect that occurs when the condition is met
			integrator.u[3] = 0
		end

		cbOn = DiscreteCallback(conditionOn,affectOn!)
		cbOff = DiscreteCallback(conditionOff,affectOff!)
		cbs = CallbackSet(cbOn,cbOff)

		for sim = 1:num_sims
			prob = DiscreteProblem([X₀; 0; 1],(0.0,fin_time),rates)
			jump_prob = JumpProblem(LV_model,prob,Direct(),save_positions=(false,false))
			sol = solve(jump_prob,SSAStepper(),tstops=protocol,saveat=0.1,callback=cbs)
			result += sum(sol[1:2,:])/((10*fin_time+1)*num_sims)
		end
		output_mean[count,:] = [durOn durOff result]
		global count += 1
	end
end

# rescale to log10
#output_mean[:,3] = log10.(output_mean[:,3])

#R code to generate heatmaps for the time average of the bacterial load
myalpha = α
mygamma = Int(log10(γ))
mykappa = κ
mymu = Int(log10(μ))
mymuprime = Int(log10(μ*m))
myxi = ξ
output_mean[:,3] = output_mean[:,3]/output_alwayson
maxSwitching = maximum(output_mean[:,3])
output_alwayson = 1

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
using RCall
@rput output_mean myalpha mygamma mykappa mymu mymuprime myxi output_alwayson maxSwitching
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
labs(x = "Time on (hours)", y = "Time off (hours)") +
scale_y_continuous(breaks=c(8,16,24,32,40),labels=c(2,4,6,8,10),expand = c(0,0)) +
scale_x_continuous(breaks=c(8,16,24,32,40),labels=c(2,4,6,8,10),expand = c(0,0)) +
ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==10^.(mymuprime)))

save_plot(p_mean,filename=paste("/home/bmorsky/antibiotic/results/mu/",paste("mu",mymu,mymuprime,sep="_"),".png", sep=""),base_height = 4,base_width = 4)
save(output_mean,file=paste("/home/bmorsky/antibiotic/results/mu/",paste("mu",mymu,mymuprime,sep="_"),".Rda", sep=""))
"""
#ggtitle(bquote(gamma==10^{-5}))
#ggtitle(bquote(gamma==2*"\u00D7"*10^{-5}))
#ggtitle(bquote(gamma==0.5*"\u00D7"*10^{-5}))
#ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==mu))
#ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==10^.(mymuprime)))
#ggtitle(bquote(kappa==.(mykappa)))
#ggtitle(bquote(mu==.(mymu)~"and"~mu^{"'"}==10^.(myM)))
#ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==10^.(mymuprime)))
#ggtitle(bquote(xi==.(myxi)))
