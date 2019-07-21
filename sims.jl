# Simulations
using Plots, Distributions

# Parameters/variables
T = 480; # Number of turns
period = 3;
num_per = T/(2*period);
t = collect(0:T);
num_sims = 100;
output1 = zeros(2*(T+1)*num_sims,4);
total_output1 = zeros((T+1)*num_sims,3);
output2 = zeros(2*(T+1)*num_sims,4);
total_output2 = zeros((T+1)*num_sims,3);
Pop = 900000;
K = 1000000.0;
m = 0.00015; # Mutation rate
b = 0.1;#0.3; # plasmid transfer rate
d = 0.85;#0.9; # Death rate
w_r = 1;#1.04; # Growth rate of the resistant type in the antibiotic regime
w = 1.01;#1.05; # Growth rate in the antibiotic-free regime
s = 0.02;#0.02; # Selective advantage of the susceptible type in the antibiotic-free regime
work1 = 0;
work2 = 0;
I = [1.0, 0.5];
I_r = [2.0, 1.0];
Ia = [1.0, 2.0];
Ia_r = [0.5, 1.0];

function antibioticOn(period,output)
	N = zeros(period+1,2);
	N[1,1] = output[1];
	N[1,2] = output[2];
	for t = 2:period+1 #for(i in 2:gens+1)	{  #simulate W-F model starting with init.num mutants
		N[t,1] = rand(Poisson(max(d*N[t-1,1]*(1 - dot(N[t,:],Ia)/K),0)));
		N[t,2] = rand(Poisson(max(w_r*N[t-1,2]*(1 - dot(N[t,:],Ia_r)/K),0))); #generate the next j as a single binomial random variable with parameters N and p.star
		if N[t,2] > 1000000
			N[t,2] = 1000000;
		end
		if N[t,1] > 0
			mutants = rand(Binomial(N[t,1], m + (1-m)*b*N[t,1]*N[t,2]/(N[t,1]+N[t,2])^2));
			N[t,1] = N[t,1] - mutants;
			N[t,2] = N[t,2] + mutants;
		end
	end
	return N[2:end,:]
end

function antibioticOff(period,output)
	N = zeros(period+1,2);
	N[1,1] = output[1];
	N[1,2] = output[2];
	for t = 2:period+1 #for(i in 2:gens+1)	{  #simulate W-F model starting with init.num mutants
		N[t,1] = rand(Poisson(max((w+s)*N[t-1,1]*(1 - dot(N[t,:],I)/K),0)));
		N[t,2] = rand(Poisson(max(w*N[t-1,2]*(1 - dot(N[t,:],I_r)/K),0))); #generate the next j as a single binomial random variable with parameters N and p.star
		if N[t,2] > 1000000
			N[t,2] = 1000000;
		end
		if N[t,2] > 0
			mutants = rand(Binomial(N[t,2], m));
			N[t,1] = N[t,1] + mutants;
			N[t,2] = N[t,2] - mutants;
		end
	end
	return N[2:end,:]
end
	n1 = zeros(1,2);
	n2 = zeros(1,2);
#for k = 2:T #this loops over the number of replicate runs

for sim = 1:num_sims
	n1 = zeros(1,2);
	n2 = zeros(1,2);
	n1[1,1] = Pop; # Initial population
	n1[1,2] = 0;
	n2[1,1] = Pop; # Initial population
	n2[1,2] = 0;
	for p = 1:num_per
		#if sum(output[end,:]) == 0
		#	break;
		#end
		n1 = vcat(n1, antibioticOn(3,n1[end,:]));
		n1 = vcat(n1, antibioticOff(3,n1[end,:]));
		n2 = vcat(n2, antibioticOn(6,n2[end,:]));
	end
	work1 += sum(n1[end,:])/1000;
	work2 += sum(n2[end,:])/1000;
	output1[(2*T+2)*(sim-1)+1:(2*T+2)*sim,:] = [vcat(t,t) reshape(n1,(2*size(n1,1),1)) vcat((2*sim-1)*ones(size(n1,1),1),2*sim*ones(size(n1,1),1)) vcat(zeros(size(n1,1),1),ones(size(n1,1),1))];
	total_output1[(T+1)*(sim-1)+1:(T+1)*sim,:] = [t sum(n1,2) sim*ones(size(n1,1),1)];
	total_output1[total_output1[:,2].>1000,2] .= 1000;
	output2[(2*T+2)*(sim-1)+1:(2*T+2)*sim,:] = [vcat(t,t) reshape(n2,(2*size(n2,1),1)) vcat((2*sim-1)*ones(size(n2,1),1),2*sim*ones(size(n2,1),1)) vcat(zeros(size(n2,1),1),ones(size(n2,1),1))];
	total_output2[(T+1)*(sim-1)+1:(T+1)*sim,:] = [t sum(n2,2) sim*ones(size(n2,1),1)];
	total_output2[total_output2[:,2].>1000,2] .= 1000;
end

using RCall
@rput output1 total_output1 output2 total_output2;
R"""
library(ggplot2)
library(cowplot)
library(extrafont)
loadfonts()
library("reshape2")

output1 <- as.data.frame(output1)
total_output1 <- as.data.frame(total_output1)
output2 <- as.data.frame(output2)
total_output2 <- as.data.frame(total_output2)

ts <- ggplot() + theme(legend.position="none") + labs(x = "time", y = "Bacterial load") + scale_x_continuous(expand = c(0, 0),limits = c(0,505)) + scale_y_continuous(expand = c(0, 0),limits = c(0,1000000))
output1 <- ts + geom_line(alpha=0.25,data=output1,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
total_output1 <- ts + geom_line(alpha=0.25,data=total_output1,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("black"),guide=FALSE)
output2 <- ts + geom_line(alpha=0.25,data=output2,aes(x=V1,y=V2,group=factor(V3),color=factor(V4))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("blue","red"),guide=FALSE) + aes(color = V3)
total_output2 <- ts + geom_line(alpha=0.25,data=total_output2,aes(x=V1,y=V2,group=factor(V3))) + ggtitle("Fixed length pulses") + scale_color_manual(values=c("black"),guide=FALSE)


plot_out <- plot_grid(output1,output2,labels=letters[1:2],label_fontfamily="LM Roman 10",ncol=2)
save_plot(plot_out,filename="~/Documents/Notre Dame/ND paper 2/Code/ts_fixed_pulsed_noK.png",base_height = 5,base_width = 10)
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
