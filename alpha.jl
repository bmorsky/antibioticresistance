using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, OrdinaryDiffEq, Statistics, ModelingToolkit

# generations where 1 generation = 1/4 hr
# Parameters
ξ = 0 # stochasticity ξ
b₁ = ξ+0.35 #0.3# susceptible birth rate
b₂ = ξ+0.3 #0.25# resistant birth rate
d = ξ+0.1 #0.05# death rate
fin_time = 1344 # length of a simulation
m = 10 # mutant stress
num_sims = 50 # number of realisations
κ = 4 # competition parameter κ
γ = 1e-5
X₀ = 1/γ
γ₁₁ = 2*γ # susceptible death rate from other susceptible
γ₁₂ = γ/κ # susceptible death rate from resistant
γ₂₁ = γ*κ # resistant death rate from susceptible
γ₂₂ = 2*γ # resistant death rate from other resistant
μ = 1e-5 # mutation rate μ

# Output
output = zeros(2*num_sims*61,3)

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

# Switching protocols
# Controller conditions
durOn = 8
durOff = 8
protocolOn = collect(0:durOn+durOff:fin_time)
protocolOff = collect(durOn:durOn+durOff:fin_time)
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

# vary alpha
count = 1
for α = 0.3:0.005:0.6
	rates = (b₁, b₂, d, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α, μ, m)

	# Always on
	for sim = 1:num_sims
		prob = DiscreteProblem([X₀; 0; 1],(0.0,fin_time),rates)
		jump_prob = JumpProblem(LV_model,prob,Direct(),save_positions=(false,false))
		sol = solve(jump_prob,SSAStepper(),saveat=0.1)
		output[count,:] = [α sum(sol[1:2,:])/(10*fin_time+1) 0]
		global count += 1
	end

	for sim = 1:num_sims
		prob = DiscreteProblem([X₀; 0; 1],(0.0,fin_time),rates)
		jump_prob = JumpProblem(LV_model,prob,Direct(),save_positions=(false,false))
		sol = solve(jump_prob,SSAStepper(),tstops=protocol,saveat=0.1,callback=cbs)
		output[count,:] = [α sum(sol[1:2,:])/(10*fin_time+1) 1]
		global count += 1
	end
end

#R code to generate alpha vs relative average bacterial load
using RCall
@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)
names(output) <- c("Parameter","Result","Type")

output <- ggplot(data=output,aes(x=Parameter,y=Result,group=Type)) + geom_point(size=1) +
stat_summary(fun=mean, geom="line") +
labs(x = bquote(alpha), y = bquote("Bacterial load"~(10^4))) +
scale_color_manual(values=c("red","blue"),labels = c("Constant", "Pulsed")) +
aes(color = factor(Type)) +
theme(legend.title=element_blank(),legend.position="right",legend.direction='vertical') +
scale_y_continuous(breaks=c(0,5000,10000,15000,20000),labels=c(0,0.5,1,1.5,2),expand = c(0, 0),limits=c(0,20100)) +
scale_x_continuous(breaks=c(0.3,0.4,0.5,0.6),expand = c(0,0),limits=c(0.3,0.603))

save_plot(output,filename="/home/bmorsky/antibiotic/results/alpha.png",base_height = 2.5,base_width = 8)
save(output,file="/home/bmorsky/antibiotic/results/alpha.Rda")
"""
#theme(legend.title=element_blank(),legend.position=c(.45,.3))
#ggtitle(bquote(gamma==10^{-5}))
#ggtitle(bquote(gamma==2*"\u00D7"*10^{-5}))
#ggtitle(bquote(gamma==0.5*"\u00D7"*10^{-5}))
#ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==mu))
#ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==10^.(mymuprime)))
#ggtitle(bquote(kappa==.(mykappa)))
#ggtitle(bquote(mu==.(mymu)~"and"~mu^{"'"}==10^.(myM)))
#ggtitle(bquote(mu==10^.(mymu)~"and"~paste(mu,"\u2032")==10^.(mymuprime)))
#ggtitle(bquote(xi==.(myxi)))
