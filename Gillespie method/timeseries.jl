# code to generate time series (Figure 1)
using Catalyst, DiffEqBase, DiffEqJump, DifferentialEquations, Plots

γ = 1e-5 # competition rate
κ = 4 # strength of competition
time=1344.0 # generations where 1 generation = 15 minutes
X₀ = 1/γ

# Parameters
b₁ = 0.35 # susceptible birth rate
b₂ = 0.3 # resistant birth rate
d = 0.1 # death rate
γ₁₁ = γ # susceptible death rate from other susceptible
γ₁₂ = γ/κ # susceptible death rate from resistant
γ₂₁ = γ*κ # resistant death rate from susceptible
γ₂₂ = γ # resistant death rate from other resistant
α = 0.4 # susceptible death rate from antibiotic
μ = 1e-5 # mutation rate
λ = 0.0
rates = [b₁, b₂, d, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α, μ, λ]

# Stochastic LV chemical reaction system
LV_model = @reaction_network begin
	b₁, n₁ → 2n₁ # n₁ birth
	b₂, n₂ → 2n₂ # n₂ birth
    d, n₁ → 0 # n₁ death
	d, n₂ → 0 # n₂ death
	2*γ₁₁, n₁+n₁ → n₁ # self competition of n₁
	γ₁₂, n₁+n₂ → n₂ # competition to n₁ from n₂
 	γ₂₁, n₁+n₂ → n₁ # competition to n₂ from n₁
	2*γ₂₂, n₂+n₂ → n₂ # self competition of n₂
	0.0001*A*α, n₁ → 0 # death by antibiotic
	0.0001*A*α/10, n₂ → 0 # death by antibiotic, α₂ < α₁
	(0.0001*9*A+1)*μ, n₁ → n₂ # mutation
	μ, n₂ → n₁ # mutation
	λ, A → 0
end b₁ b₂ d γ₁₁ γ₁₂ γ₂₁ γ₂₂ α μ λ

# Solve system for various protocols averaged 100 times
# Callbacks
durOn = 20#48 # 2.5hrs on
durOff = 20#48 # 2.5hrs off

# Controller conditions
protocolOn = collect(0:durOn+durOff:time)
protocolOff = collect(durOn:durOn+durOff:time)
# protocol = merge(protocolOn,protocolOff)
protocol = collect(0:1:time)
conditionOn(u,t,integrator) = t in (protocolOn)
conditionOff(u,t,integrator) = t in (protocolOff)

function affectOn!(integrator) # effect that occurs when the condition is met
	integrator.p[10] = 0
	integrator.u[3] = 10000
	reset_aggregated_jumps!(integrator)
end
function affectOff!(integrator) # effect that occurs when the condition is met
	integrator.p[10] = 0.1
	reset_aggregated_jumps!(integrator)
end
cbOn = DiscreteCallback(conditionOn,affectOn!,save_positions = (false,false))
cbOff = DiscreteCallback(conditionOff,affectOff!,save_positions = (false,false))
cbs = CallbackSet(cbOn,cbOff)

prob = DiscreteProblem([X₀; 0; 10000],(0.0,time),rates)
jump_prob = JumpProblem(LV_model,prob,Direct(),save_positions=(false,false))
sol = solve(jump_prob,SSAStepper(),saveat=0.1,tstops=protocol,callback=cbs)
output1 = [sol[3,:] sol.t[:] zeros(length(sol[1,:]),1); sol[2,:] sol.t[:] ones(length(sol[2,:]),1);] #output

prob = DiscreteProblem([X₀; 0; 10000],(0.0,time),rates)
jump_prob = JumpProblem(LV_model,prob,Direct(),save_positions=(false,false))
sol = solve(jump_prob,SSAStepper(),saveat=0.1)
output2 = [sol[1,:] sol.t[:] zeros(length(sol[1,:]),1); sol[2,:] sol.t[:] ones(length(sol[2,:]),1);] #output

#R code to generate figure
using RCall
@rput output1 output2 time
R"""
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot());

output1 <- as.data.frame(output1)
names(output1) <- c("Number","Time","Type")
output2 <- as.data.frame(output2)
names(output2) <- c("Number","Time","Type")

output1 <- ggplot() + geom_line(data=output1,aes(x=Time,y=Number,group=Type),size=1) +
ggtitle("On/off protocol") + labs(x = "Time (days)", y = bquote("Bacterial load"~(10^4))) +
scale_color_manual(values=c("blue","red"),labels = c("Wild-type", "Mutant")) +
aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,.65)) +
scale_y_continuous(breaks=c(0,10000,20000,30000,40000),labels=c(0,1,2,3,4),expand = c(0, 0),limits=c(0,40000)) +
scale_x_continuous(breaks=c(0,192,384,576,768,960,1152,1344),labels=c(0,2,4,6,8,10,12,14),expand = c(0,0))

output2 <- ggplot() + geom_line(data=output2,aes(x=Time,y=Number,group=Type),size=1) +
ggtitle("Constant application") + labs(x = "Time (days)", y = bquote("Bacterial load"~(10^4))) +
scale_color_manual(values=c("blue","red"),labels = c("Wild-type", "Mutant")) +
aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,.65)) +
scale_y_continuous(breaks=c(0,10000,20000,30000,40000),labels=c(0,1,2,3,4),expand = c(0, 0),limits=c(0,40000)) +
scale_x_continuous(breaks=c(0,192,384,576,768,960,1152,1344),labels=c(0,2,4,6,8,10,12,14),expand = c(0,0))

save_plot(output1,filename="~/Desktop/ts24hrperiod.png",base_height = 2.5,base_width = 10)
save_plot(output2,filename="~/Desktop/tsOn.png",base_height = 2.5,base_width = 10)
"""
