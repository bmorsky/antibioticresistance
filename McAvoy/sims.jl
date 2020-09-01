# Simulations for the stochastic model for a given protocol. Results plotted as
# a time series for the protocol vs a constant antibiotic application.
using Plots, Distributions, LinearAlgebra

#Parameters
d = 0.9 #death rate
I = [1.0, 1.0] #susceptible type interspecific competition parameters (antibiotic off)
I_r = [1.0, 1.0] #resistant type interspecific competition parameters (antibiotic off)
Ia = [1.0, 1.0] #susceptible type interspecific competition parameters (antibiotic on)
Ia_r = [1.0, 1.0] #resistant type interspecific competition parameters (antibiotic on)
init_pop = 900 #initial population size
K = 1000.0 #carrying capacity
m = 0.0002 #mutation rate
num_sims = 100 #number of simulations per parameter combination
w_r = 1.02 #growth rate of the resistant type in the antibiotic regime
w = 1.02 #growth rate in the antibiotic-free regime
s = 0.02 #elective advantage of the susceptible type in the antibiotic-free regime

#Protocol parameters
dur_off = 3 #duration the antibiotic is off
dur_on = 3 #duration the antibiotic is on
num_periods = 80 #number of periods
num_turns = num_periods*(dur_on+dur_off) #number of turns

#Output
t = collect(0:num_turns)
output1 = zeros(2*(num_turns+1)*num_sims,4)
total_output1 = zeros((num_turns+1)*num_sims,3)
output2 = zeros(2*(num_turns+1)*num_sims,4)
total_output2 = zeros((num_turns+1)*num_sims,3)

#Variables to track
work1 = 0
work2 = 0
n1 = zeros(1,2) #number of susceptible and resistant bacteria for the protocol
n2 = zeros(1,2) #number of susceptible and resistant bacteria for the constant antibiotic application case

#Dynamics for when the antibiotic is on
function antibioticOn(period,output)
	N = zeros(period+1,2)
	N[1,1] = output[1]
	N[1,2] = output[2]
	for t = 2:period+1
		N[t,1] = rand(Poisson(max(d*N[t-1,1]*(1 - dot(N[t,:],Ia)/K),0)))
		N[t,2] = rand(Poisson(max(w_r*N[t-1,2]*(1 - dot(N[t,:],Ia_r)/K),0)))
		if N[t,2] > 1000000
			N[t,2] = 1000000
		end
		if N[t,1] > 0
			mutants = rand(Binomial(N[t,1], m))
			N[t,1] = N[t,1] - mutants
			N[t,2] = N[t,2] + mutants
		end
	end
	return N[2:end,:]
end

#Dynamics for when the antibiotic is off
function antibioticOff(period,output)
	N = zeros(period+1,2)
	N[1,1] = output[1]
	N[1,2] = output[2]
	for t = 2:period+1
		N[t,1] = rand(Poisson(max((w+s)*N[t-1,1]*(1 - dot(N[t,:],I)/K),0)))
		N[t,2] = rand(Poisson(max(w*N[t-1,2]*(1 - dot(N[t,:],I_r)/K),0)))
		if N[t,2] > 1000000
			N[t,2] = 1000000
		end
		if N[t,2] > 0
			mutants = rand(Binomial(N[t,2], m))
			N[t,1] = N[t,1] + mutants
			N[t,2] = N[t,2] - mutants
		end
	end
	return N[2:end,:]
end

#Runs simulations
for sim = 1:num_sims
	#Initialize population
	n1 = zeros(1,2)
	n2 = zeros(1,2)
	n1[1,1] = init_pop #initial number of susceptible type for the protocol
	n1[1,2] = 0 #initial number of resistant type for the protocol
	n2[1,1] = init_pop #initial number of susceptible type for the constant antibiotic application case
	n2[1,2] = 0 #initial number of resistant type for the constant antibiotic application case
	for p = 1:num_periods
		#if sum(output[end,:]) == 0
		#	break;
		#end
		n1 = vcat(n1, antibioticOn(dur_on,n1[end,:])) #apply antibiotic for dur_on turns
		n1 = vcat(n1, antibioticOff(dur_off,n1[end,:])) #no antibiotic for dur_off turns
		n2 = vcat(n2, antibioticOn(dur_on+dur_off,n2[end,:])) #apply antibiotic for dur_on+dur_off turns
	end
	global work1 += sum(n1[end,:])/1000
	global work2 += sum(n2[end,:])/1000
	output1[(2*num_turns+2)*(sim-1)+1:(2*num_turns+2)*sim,:] = [vcat(t,t) reshape(n1,(2*size(n1,1),1)) vcat((2*sim-1)*ones(size(n1,1),1),2*sim*ones(size(n1,1),1)) vcat(zeros(size(n1,1),1),ones(size(n1,1),1))]
	total_output1[(num_turns+1)*(sim-1)+1:(num_turns+1)*sim,:] = [t sum(n1,dims=2) sim*ones(size(n1,1),1)]
	total_output1[total_output1[:,2].>1000,2] .= 1000
	output2[(2*num_turns+2)*(sim-1)+1:(2*num_turns+2)*sim,:] = [vcat(t,t) reshape(n2,(2*size(n2,1),1)) vcat((2*sim-1)*ones(size(n2,1),1),2*sim*ones(size(n2,1),1)) vcat(zeros(size(n2,1),1),ones(size(n2,1),1))]
	total_output2[(num_turns+1)*(sim-1)+1:(num_turns+1)*sim,:] = [t sum(n2,dims=2) sim*ones(size(n2,1),1)]
	total_output2[total_output2[:,2].>1000,2] .= 1000
end

#R code to generate figure
using RCall
@rput output1 total_output1 output2 total_output2;
R"""
library(ggplot2)
library(cowplot)

output1 <- as.data.frame(output1)
total_output1 <- as.data.frame(total_output1)
output2 <- as.data.frame(output2)
total_output2 <- as.data.frame(total_output2)

ts <- ggplot() + theme(legend.position="none") + labs(x = "time", y = "Bacterial load") + scale_x_continuous(expand = c(0, 0),limits = c(0,505)) + scale_y_continuous(expand = c(0, 0),limits = c(0,1000))
output1 <- ts + geom_line(alpha=0.25,data=output1,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Protocol") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
total_output1 <- ts + geom_line(alpha=0.25,data=total_output1,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Protocol") + scale_color_manual(values=c("black"),guide=FALSE)
output2 <- ts + geom_line(alpha=0.25,data=output2,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Constant application") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
total_output2 <- ts + geom_line(alpha=0.25,data=total_output2,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Constant application") + scale_color_manual(values=c("black"),guide=FALSE)

plot_out <- plot_grid(output1,output2,labels=letters[1:2],ncol=2)
save_plot(plot_out,filename="~/Documents/Notre Dame/ND paper 2/Code/antibioticresistance/McAvoy/protVconst.png",base_height = 5,base_width = 10)
"""

#p = N[t-1,1]*(1+s)*(1-m)/(N[t-1,1]*(1+s)+N[t-1,2]);  #compute the post-selection expected frequency
# This simualates the Wright-Fisher Model of selection and random genetic drift
# for an asexual organism with two alleles and no mutation
# By R. Gomulkiewicz (5 September 2013)
# Updated 8 Sep 2015

#s <- 0.01 		# selection coefficient of the mutant typeS
#init.num <- 1	# initial number of the mutant type
#N <-10		# population size N
#reps <-	15		# number of replicate simulation runs
#gens <- 40 		# generations of evolution to simulate for each replicate run

#create a matrix with gens+1 rows  and reps columns whose every entry is the initial number of mutants
#below,we will replace entries 2 through gens+1 of each column with simulated numbers of mutants in each replicate run
#j=matrix(init.num,gens+1,reps)

#loops that execute the replicate simulations
#uses simplified form of the post-selection expected frequency pstar = j(1+s)/[j(1+s)+(N-j)] = j(1+s)/(js + N)
#for(k in 1:reps){ #this loops over the number of replicate runs
#for(i in 2:gens+1)	{  #simulate W-F model starting with init.num mutants
#	p.star <- j[i-1,k]*(1+s)/(j[i-1,k]*s+N)  #compute the post-selection expected frequency, given j
#	j[i,k]=rbinom(1,N,p.star)  #generate the next j as a single binomial random variable with parameters N and p.star
#	}
#}

#Show the results as mutant frequencies (i.e., plot each mutant count divided by 2N)
#Display all the replicate runs on a single plot
#matplot(0:gens,j/N,type="l",ylab="freq(A)",xlab="generation")