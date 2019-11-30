# Optimal control for the mean field model. Finds the optimal protocol for an
# initial condition of 900 susceptible bacteria by solving the state and adjoint
# equations. However, it does not prevent the emergence and establishment of the
# resistant type. Plots the time series of susceptible and resistance bacteria
# and the protocol.
using Plots, Distributions, DifferentialEquations, RecursiveArrayTools, Sundials

#Parameters
β₁ = 0.02 # susceptible birth rate
δ₁ = 0.014 # susceptible death rate
β₂ = 0.01 # resistant birth rate
δ₂ = 0.005 # resistant death rate
γ₁₁ = 1e-6 # susceptible death rate from other susceptible
γ₁₂ = 1e-6 # susceptible death rate from resistant
γ₂₁ = 1e-6 # resistant death rate from susceptible
γ₂₂ = 1e-6 # resistant death rate from other resistant
α₁ = 0.008 # susceptible death rate from antibiotic
α₂ = 0.0015 # resistant death rate from antibiotic
μ = 0.10e-5 # mutation rate
μₐ = 10*μ # mutation rate cause by antibiotic
η = 0.0 #4e-3 # pasmid transfer rate
init_pop = 1e4 #initial population size
rates = (β₁, δ₁, β₂, δ₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ, μₐ, η)

#System of ODEs: u[1] susceptible type, u[2] resistant type, u[3] adjoit equation for u[1], adjoint equation for u[2], control u[5].
function antiOn(du,u,p,t)
	du[1] = (β₁-δ₁-μ)*u[1] - γ₁₁*u[1]^2 - (γ₁₂+η)*u[1]*u[2] + μ*u[2] - (α₁+μₐ)*u[1]
	du[2] = (β₂-δ₂-μ)*u[2] - γ₂₂*u[2]^2 + (η-γ₂₁)*u[1]*u[2] + μ*u[1] + (μₐ*u[1] - α₂*u[2])
end #β₁ δ₁ β₂ δ₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ η
function antiOff(du,u,p,t)
	du[1] = (β₁-δ₁-μ)*u[1] - γ₁₁*u[1]^2 - (γ₁₂+η)*u[1]*u[2] + μ*u[2]
	du[2] = (β₂-δ₂-μ)*u[2] - γ₂₂*u[2]^2 + (η-γ₂₁)*u[1]*u[2] + μ*u[1]
end #β₁ δ₁ β₂ δ₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ μₐ η


# Run multiple realizations. Solve system num_sims times and plot the results.
#for init = 0.0:200.0:1000.0
	u0 = [100.0,0.0] #initial conditions: init_pop susceptible, 0 resistance, adjoint equations 1 and 1 for control minimizing the total number of bacteria at the final time, and the antibiotic on i.e. a = 1
 	tspan = (0,10000.0) #time span
 	prob = ODEProblem(antiOff,u0,tspan,rates) #the problem to solve
  	sol = solve(prob,Tsit5())
#  	plot!(sol,vars = (1,2))
# end
# plot!(legend=false,vars = (1,2))
sol[end]
 #System of ODEs: u[1] susceptible type, u[2] resistant type, u[3] adjoit equation for u[1], adjoint equation for u[2], control u[5].

 # # Run multiple realizations. Solve system num_sims times and plot the results.
 # for init = 0.0:100.0:1000.0
 # 	u0 = [init,1000.0-init] #initial conditions: init_pop susceptible, 0 resistance, adjoint equations 1 and 1 for control minimizing the total number of bacteria at the final time, and the antibiotic on i.e. a = 1
 # 	tspan = (0,100000.0) #time span
 # 	prob = ODEProblem(antiOff,u0,tspan,rates) #the problem to solve
 #  	sol = solve(prob,Tsit5())
 #  	plot!(sol,vars = (1,2))
 # end
 # plot!(legend=false,vars = (1,2))

#for plotting in R
output = [sol[1,:] sol.t[:] zeros(length(sol[1,:]),1); sol[2,:] sol.t[:] ones(length(sol[2,:]),1)] #output
#protocol = [sol[5,:] sol.t[:]] #the protocol
#R code to generate figure
using RCall
@rput output;
R"""
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

output <- as.data.frame(output)
names(output) <- c("Number","Time","Type")

ts <- ggplot() + scale_x_continuous(expand = c(0, 0))

output <- ts + geom_line(data=output,aes(x=Time,y=Number,group=Type),size=1) + ggtitle("Time series for the optimal control")+ scale_color_manual(values=c("blue","red"),labels = c("Susceptible", "Resistant")) + aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,.85)) + labs(x = "Time", y = "Bacterial load") + scale_y_continuous(expand = c(0, 0))

plot_out <- plot_grid(output,labels=letters[1:2],ncol=2)
save_plot(plot_out,filename="~/Documents/Notre Dame/ND paper 2/Code/phase.png",base_height = 5,base_width = 10)
"""

################ OLD CODE BELOW #####################

# mutable struct SimType{T} <: DEDataVector{T}
#     x::Array{T,1}
#     f1::T
# end
#
# #f = @ode_def LotkaVolterra begin
# #  dx = -0.5*u + 2*u.f1
# #  dy = -0.5*u
# #end
#
# function f(du,u,p,t)
#     du[1] = -0.5*u[1] + u.f1
#     du[2] = -0.5*u[2]
# end
# const tstop1 = [5.]
# const tstop2 = [8.]
#
#
# function condition(u,t,integrator)
#   t in tstop1
# end
#
# function condition2(u,t,integrator)
#   t in tstop2
# end
# function affect!(integrator)
#   for c in full_cache(integrator)
#     c.f1 = 1
#   end
# end
#
# function affect2!(integrator)
#   for c in full_cache(integrator)
#     c.f1 = 0
#   end
# end
# save_positions = (true,true)
#
# cb = DiscreteCallback(condition, affect!, save_positions=save_positions)
#
# save_positions = (false,true)
#
# cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)
#
# cbs = CallbackSet(cb,cb2)
# u0 = SimType([10.0;10.0], 0.0)
# prob = ODEProblem(f,u0,(0.0,10.0))
# const tstop = [5.;8.]
# sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
# plot(sol)
#
#
#
# antiOn = @ode_def LotkaVolterra begin
#   dx = d*x*(1 + (x+y)/K) - m*x - x*(1-m)*b*y*(x+y)
#   dy = wprime*x*(1 - (x+y)/K) + m*x + x*(1-m)*b*y*(x+y)
# end d K m b wprime
#
# antiOff = @ode_def LotkaVolterra begin
#   dx = (w+s)*(1 + (x+y)/K)
#   dy = w*(1 - (x+y)/K)
# end w s K
#
# u0 = [50.0,0.0]
# pOn = (-0.9,1000.0,0.001,0.01,1.2)
# pOff = (-0.9,1000.0,0.001,0.01,1.2)
# tspan = (0.0,0.25)
#
# out = zeros(1,2)
#
# for t=1:10
# 	prob = ODEProblem(antiOn,u0,tspan,pOn)
# 	sol = solve(prob)
# 	out = vcat(out,sol)
# 	u0 = sol[end]
# 	tspan = (0.25*t, 0.25*(t+1))
# 	prob = ODEProblem(antiOff,u0,tspan,pOff)
# 	sol = solve(prob)
# 	u0 = sol[end]
# 	tspan = (0.25*(t+1), 0.25*(t+2))
# 	out = vcat(out,sol)
# end
#
# plot(sol)
#
# # Parameters/variables
# T = 480; # Number of turns
# period = 3;
# num_per = T/(2*period);
# t = collect(0:T);
# num_sims = 100;
# output1 = zeros(2*(T+1)*num_sims,4);
# total_output1 = zeros((T+1)*num_sims,3);
# output2 = zeros(2*(T+1)*num_sims,4);
# total_output2 = zeros((T+1)*num_sims,3);
# Pop = 900000;
# K = 1000000.0;
# m = 0.00015; # Mutation rate
# b = 0.1;#0.3; # plasmid transfer rate
# d = 0.85;#0.9; # Death rate
# w_r = 1;#1.04; # Growth rate of the resistant type in the antibiotic regime
# w = 1.01;#1.05; # Growth rate in the antibiotic-free regime
# s = 0.02;#0.02; # Selective advantage of the susceptible type in the antibiotic-free regime
# work1 = 0;
# work2 = 0;
# I = [1.0, 0.5];
# I_r = [2.0, 1.0];
# Ia = [1.0, 2.0];
# Ia_r = [0.5, 1.0];
#
# function antibioticOn(period,output)
# 	N = zeros(period+1,2);
# 	N[1,1] = output[1];
# 	N[1,2] = output[2];
# 	for t = 2:period+1 #for(i in 2:gens+1)	{  #simulate W-F model starting with init.num mutants
# 		N[t,1] = rand(Poisson(max(d*N[t-1,1]*(1 - dot(N[t,:],Ia)/K),0)));
# 		N[t,2] = rand(Poisson(max(w_r*N[t-1,2]*(1 - dot(N[t,:],Ia_r)/K),0))); #generate the next j as a single binomial random variable with parameters N and p.star
# 		if N[t,2] > 1000000
# 			N[t,2] = 1000000;
# 		end
# 		if N[t,1] > 0
# 			mutants = rand(Binomial(N[t,1], m + (1-m)*b*N[t,1]*N[t,2]/(N[t,1]+N[t,2])^2));
# 			N[t,1] = N[t,1] - mutants;
# 			N[t,2] = N[t,2] + mutants;
# 		end
# 	end
# 	return N[2:end,:]
# end
#
# function antibioticOff(period,output)
# 	N = zeros(period+1,2);
# 	N[1,1] = output[1];
# 	N[1,2] = output[2];
# 	for t = 2:period+1 #for(i in 2:gens+1)	{  #simulate W-F model starting with init.num mutants
# 		N[t,1] = rand(Poisson(max((w+s)*N[t-1,1]*(1 - dot(N[t,:],I)/K),0)));
# 		N[t,2] = rand(Poisson(max(w*N[t-1,2]*(1 - dot(N[t,:],I_r)/K),0))); #generate the next j as a single binomial random variable with parameters N and p.star
# 		if N[t,2] > 1000000
# 			N[t,2] = 1000000;
# 		end
# 		if N[t,2] > 0
# 			mutants = rand(Binomial(N[t,2], m));
# 			N[t,1] = N[t,1] + mutants;
# 			N[t,2] = N[t,2] - mutants;
# 		end
# 	end
# 	return N[2:end,:]
# end
# 	n1 = zeros(1,2);
# 	n2 = zeros(1,2);
# #for k = 2:T #this loops over the number of replicate runs
#
# for sim = 1:num_sims
# 	n1 = zeros(1,2);
# 	n2 = zeros(1,2);
# 	n1[1,1] = Pop; # Initial population
# 	n1[1,2] = 0;
# 	n2[1,1] = Pop; # Initial population
# 	n2[1,2] = 0;
# 	for p = 1:num_per
# 		#if sum(output[end,:]) == 0
# 		#	break;
# 		#end
# 		n1 = vcat(n1, antibioticOn(3,n1[end,:]));
# 		n1 = vcat(n1, antibioticOff(3,n1[end,:]));
# 		n2 = vcat(n2, antibioticOn(6,n2[end,:]));
# 	end
# 	work1 += sum(n1[end,:])/1000;
# 	work2 += sum(n2[end,:])/1000;
# 	output1[(2*T+2)*(sim-1)+1:(2*T+2)*sim,:] = [vcat(t,t) reshape(n1,(2*size(n1,1),1)) vcat((2*sim-1)*ones(size(n1,1),1),2*sim*ones(size(n1,1),1)) vcat(zeros(size(n1,1),1),ones(size(n1,1),1))];
# 	total_output1[(T+1)*(sim-1)+1:(T+1)*sim,:] = [t sum(n1,2) sim*ones(size(n1,1),1)];
# 	total_output1[total_output1[:,2].>1000,2] .= 1000;
# 	output2[(2*T+2)*(sim-1)+1:(2*T+2)*sim,:] = [vcat(t,t) reshape(n2,(2*size(n2,1),1)) vcat((2*sim-1)*ones(size(n2,1),1),2*sim*ones(size(n2,1),1)) vcat(zeros(size(n2,1),1),ones(size(n2,1),1))];
# 	total_output2[(T+1)*(sim-1)+1:(T+1)*sim,:] = [t sum(n2,2) sim*ones(size(n2,1),1)];
# 	total_output2[total_output2[:,2].>1000,2] .= 1000;
# end
#
# using RCall
# @rput output1 total_output1 output2 total_output2;
# R"""
# library(ggplot2)
# library(cowplot)
# library(extrafont)
# loadfonts()
# library("reshape2")
#
# output1 <- as.data.frame(output1)
# total_output1 <- as.data.frame(total_output1)
# output2 <- as.data.frame(output2)
# total_output2 <- as.data.frame(total_output2)
#
# ts <- ggplot() + theme(legend.position="none") + labs(x = "time", y = "Bacterial load") + scale_x_continuous(expand = c(0, 0),limits = c(0,505)) + scale_y_continuous(expand = c(0, 0),limits = c(0,1000000))
# output1 <- ts + geom_line(alpha=0.25,data=output1,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
# total_output1 <- ts + geom_line(alpha=0.25,data=total_output1,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("black"),guide=FALSE)
# output2 <- ts + geom_line(alpha=0.25,data=output2,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
# total_output2 <- ts + geom_line(alpha=0.25,data=total_output2,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("black"),guide=FALSE)
#
#
# plot_out <- plot_grid(output1,output2,labels=letters[1:2],label_fontfamily="LM Roman 10",ncol=2)
# save_plot(plot_out,filename="~/Documents/Notre Dame/ND paper 2/Code/ts_fixed_pulsed_noK.png",base_height = 5,base_width = 10)
# """
#
# 		#p = N[t-1,1]*(1+s)*(1-m)/(N[t-1,1]*(1+s)+N[t-1,2]);  #compute the post-selection expected frequency
# # This simualates the Wright-Fisher Model of selection and random genetic drift
# # for an asexual organism with two alleles and no mutation
# # By R. Gomulkiewicz (5 September 2013)
# # Updated 8 Sep 2015
#
# #s <- 0.01 		# selection coefficient of the mutant typeS
# #init.num <- 1	# initial number of the mutant type
# #N <-10		# population size N
# #reps <-	15		# number of replicate simulation runs
# #gens <- 40 		# generations of evolution to simulate for each replicate run
#
# #create a matrix with gens+1 rows  and reps columns whose every entry is the initial number of mutants
# #below,we will replace entries 2 through gens+1 of each column with simulated numbers of mutants in each replicate run
# #j=matrix(init.num,gens+1,reps)
#
# #loops that execute the replicate simulations
# #uses simplified form of the post-selection expected frequency pstar = j(1+s)/[j(1+s)+(N-j)] = j(1+s)/(js + N)
# #for(k in 1:reps){ #this loops over the number of replicate runs
# #for(i in 2:gens+1)	{  #simulate W-F model starting with init.num mutants
# #	p.star <- j[i-1,k]*(1+s)/(j[i-1,k]*s+N)  #compute the post-selection expected frequency, given j
# #	j[i,k]=rbinom(1,N,p.star)  #generate the next j as a single binomial random variable with parameters N and p.star
# #	}
# #}
#
# #Show the results as mutant frequencies (i.e., plot each mutant count divided by 2N)
# #Display all the replicate runs on a single plot
# #matplot(0:gens,j/N,type="l",ylab="freq(A)",xlab="generation")
