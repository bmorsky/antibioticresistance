# Optimal control for the mean field model. Finds the optimal protocol for an
# initial condition of 900 susceptible bacteria by solving the state and adjoint
# equations. However, it does not prevent the emergence and establishment of the
# resistant type. Plots the time series of susceptible and resistance bacteria
# and the protocol.
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
ζ = 1
η = 0.000001
pOff = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, μ, η)
pOn = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, μₐ, ζ, η)
p = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, μₐ, ζ, η)
init_pop = 500.0
num_sims = 1 #number of simulations per parameter combination
# Protocol parameters
durOn = 2.0 #duration the antibiotic is off
durOff = 2.0 #duration the antibiotic is on
num_periods = 10 #number of periods
num_turns = num_periods*(durOn+durOff) #number of turns

#Variables to track
n1 = zeros(1,2) #number of susceptible and resistant bacteria for the protocol
n2 = zeros(1,2) #number of susceptible and resistant bacteria for the constant antibiotic application case
output1 = zeros(2*(num_turns+1)*num_sims,4)

# Antibiotic off
LV_antiOff = @reaction_network begin
    b₁, n₁ → 2n₁
    d₁, n₁ → 0
	b₂, n₂ → 2n₂
	d₂, n₂ → 0
	γ₁₁, n₁ + n₁ → 0
	γ₁₂, n₁ + n₂ → 0
	γ₂₁, n₂ + n₁ → 0
	γ₂₂, n₂ + n₂ → 0
	μ, n₂ → n₁
	μ, n₂ → n₁
	η, n₁ + n₂ → 2n₂
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ μ η
# Antibiotic on
LV_antiOn = @reaction_network begin
	b₁, n₁ → 2n₁
    d₁, n₁ → 0
	b₂, n₂ → 2n₂
	d₂, n₂ → 0
	γ₁₁, n₁ + n₁ → 0
	γ₁₂, n₁ + n₂ → 0
	γ₂₁, n₂ + n₁ → 0
	γ₂₂, n₂ + n₂ → 0
	α₁, a + n₁ → a
	α₂, a + n₂ → a
	μ, n₁ → n₂
	μ, n₂ → n₁
	μₐ, a + n₁ → a + n₂
	#ζ, a → 0
	η, n₁ + n₂ → 2n₂
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ ζ η
# General reactions
LV_general = @reaction_network begin
    b₁, n₁ → 2n₁
    d₁, n₁ → 0
	b₂, n₂ → 2n₂
	d₂, n₂ → 0
	γ₁₁, n₁ + n₁ → 0
	γ₁₂, n₁ + n₂ → 0
	γ₂₁, n₂ + n₁ → 0
	γ₂₂, n₂ + n₂ → 0
	α₁, a + n₁ → a
	α₂, a + n₂ → a
	μ, n₁ → n₂
	μ, n₂ → n₁
	μₐ, a + n₁ → a + n₂
	ζ, a → 0
	η, n₁ + n₂ → 2n₂
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ ζ η

# Runs simulations
# Initialize population
global n1 = zeros(2,1)
global n2 = zeros(2,1)
n1[1,1] = init_pop #initial number of susceptible type for the protocol
n1[2,1] = 0 #initial number of resistant type for the protocol
n2[1,1] = init_pop #initial number of susceptible type for the constant antibiotic application case
n2[2,1] = 0 #initial number of resistant type for the constant antibiotic application case
time1 = [0]
function myProtocol(n1,time1)
	for p = 1:num_periods
		prob = DiscreteProblem([n1[1,end]; n1[2,end]; 2.0],(0.0,durOn),pOn)
		jump_prob = JumpProblem(prob,Direct(),LV_antiOn)
		sol = solve(jump_prob,FunctionMap())
		n1 = hcat(n1, sol[1:2,:])
		time1 = vcat(time1, sol.t[:]+time1[end]*ones(length(sol.t[:])))

		prob = DiscreteProblem([n1[1,end]; n1[2,end]],(0.0,durOff),pOff)
		jump_prob = JumpProblem(prob,Direct(),LV_antiOff)
		sol = solve(jump_prob,FunctionMap())
		n1 = hcat(n1, sol[:,:])
		time1 = vcat(time1, sol.t[:]+time1[end]*ones(length(sol.t[:])))
	end
	return n1, time1
end
n1,time1 = myProtocol(n1,time1)
prob = DiscreteProblem([init_pop; 0.0; 2.0],(0.0,num_turns),pOn)
jump_prob = JumpProblem(prob,Direct(),LV_antiOn)
sol = solve(jump_prob,FunctionMap())
plot(time1,[n1[1,:],n1[2,:]])
plot(sol.t[:],[sol[1,:],sol[2,:]])

# #R code to generate figure
# using RCall
# @rput output1 total_output1 output2 total_output2;
# R"""
# library(ggplot2)
# library(cowplot)
#
# output1 <- as.data.frame(output1)
# total_output1 <- as.data.frame(total_output1)
# output2 <- as.data.frame(output2)
# total_output2 <- as.data.frame(total_output2)
#
# ts <- ggplot() + theme(legend.position="none") + labs(x = "time", y = "Bacterial load") + scale_x_continuous(expand = c(0, 0),limits = c(0,505)) + scale_y_continuous(expand = c(0, 0),limits = c(0,1000))
# output1 <- ts + geom_line(alpha=0.25,data=output1,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Protocol") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
# total_output1 <- ts + geom_line(alpha=0.25,data=total_output1,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Protocol") + scale_color_manual(values=c("black"),guide=FALSE)
# output2 <- ts + geom_line(alpha=0.25,data=output2,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Constant application") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
# total_output2 <- ts + geom_line(alpha=0.25,data=total_output2,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Constant application") + scale_color_manual(values=c("black"),guide=FALSE)
#
# plot_out <- plot_grid(output1,output2,labels=letters[1:2],ncol=2)
# save_plot(plot_out,filename="~/Documents/Notre Dame/ND paper 2/Code/ts_protocol_vs_constant_app_test.png",base_height = 5,base_width = 10)
# """

# solve
# prob = DiscreteProblem([500,500],(0.0,200.0),p)
# jump_prob = JumpProblem(prob,Direct(),LV_antiOff)
# sol = solve(jump_prob,FunctionMap())
# plot(sol)
# prob = DiscreteProblem([1000,0,10],(0.1,1000.0),p)
# jump_prob = JumpProblem(prob,Direct(),LV_antiOn)
# sol = solve(jump_prob,FunctionMap())
# plot(sol)
# prob = DiscreteProblem([500,10,20],(0.0,100.0),p)
# jump_prob = JumpProblem(prob,Direct(),LV_general)
# sol = solve(jump_prob,FunctionMap())
# plot(sol)
