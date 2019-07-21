# Simulations
using Plots, Distributions

# Parameters/variables
num_sims = 100;
Pop = 900;
#m = 0.00008; # Mutation rate
#b = 0.3; # plasmid transfer rate
#d = 0.84; # Death rate
#w_r = 1.05; # Growth rate of the resistant type in the antibiotic regime
#w = 1; # Growth rate in the antibiotic-free regime
#s = 0.1; # Selective advantage of the susceptible type in the antibiotic-free regime

m = 0.0013; # Mutation rate
b = 0.3; # plasmid transfer rate
d = 0.854; # Death rate
w_r = 1.005; # Growth rate of the resistant type in the antibiotic regime
w = 1.055; # Growth rate in the antibiotic-free regime
s = 0.02; # Selective advantage of the susceptible type in the antibiotic-free regime

maxOn = 20;
maxOff = 20;
avg_work = zeros(maxOn*(maxOff+1),3);
std_work = zeros(maxOn*(maxOff+1),3);
avg_t = zeros(maxOn*(maxOff+1),3);
std_t = zeros(maxOn*(maxOff+1),3);
K = 1000.0;
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
		if N[t,2] > 1000
			N[t,2] = 1000;
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
		if N[t,2] > 1000
			N[t,2] = 1000;
		end
		if N[t,2] > 0
			mutants = rand(Binomial(N[t,2], m));
			N[t,1] = N[t,1] + mutants;
			N[t,2] = N[t,2] - mutants;
		end
	end
	return N[2:end,:]
end

n = zeros(1,2);
count = 1;
for perOn = 5:5:5*maxOn
	for perOff = 0:5:5*maxOff
		if perOn == 0 && perOff == 0
			break;
		end
		work = zeros(num_sims,2);
		for sim = 1:num_sims
			n = zeros(1,2);
			n[1,1] = Pop; # Initial population
			n[1,2] = 0;
			t = 0;
			while sum(n[end,:]) > 0 && t <= 100 #n[end,2] < 2000 && n[end,1] < 2000
				n = vcat(n, antibioticOn(perOn,n[end,:]));
				n = vcat(n, antibioticOff(perOff,n[end,:]));
				t += perOn + perOff;
			end
			work[sim,:] = [convert(Float64, sum(n[end,:]) == 0.0) t];
		end
		avg_work[count,:] = [perOn perOff mean(work[:,1])]
		std_work[count,:] = [perOn perOff std(work[:,1])]
		avg_t[count,:] = [perOn perOff mean(work[:,2])]
		std_t[count,:] = [perOn perOff std(work[:,2])]
		count += 1;
	end
end

using RCall
@rput avg_work std_work avg_t std_t;
R"""
library(ggplot2)
library(cowplot)
library(extrafont)
library(viridis)
loadfonts()
library("reshape2")

avg_work <- as.data.frame(avg_work)
std_work <- as.data.frame(std_work)
avg_t <- as.data.frame(avg_t)
std_t <- as.data.frame(std_t)

q <- ggplot() + theme(strip.background=element_blank()) + theme(plot.margin = unit(c(0, 0, 1.25, 0), "cm")) + coord_equal() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(text = element_text(family = "LM Roman 10"))+guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top"))

q_avg_work <- q + geom_raster(data=avg_work,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="magma",limits = c(0,1)) + labs(x = "On", y = "Off", fill = "Mean success")
q_std_work <- q + geom_raster(data=std_work,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="magma",limits = c(0,1)) + labs(x = "On", y = "Off", fill = "St dev success")
q_avg_t <- q + geom_raster(data=avg_t,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="magma",limits = c(0,150)) + labs(x = "On", y = "Off", fill = "Mean time")
q_std_t <- q + geom_raster(data=std_t,aes(x=V1,y=V2,fill=V3)) + scale_fill_viridis(option="magma",limits = c(0,150)) + labs(x = "On", y = "Off", fill = "St dev time")

save_plot(plot_grid(q_avg_work,q_std_work,q_avg_t,q_std_t,label_fontfamily="LM Roman 10",ncol=2),filename="~/Documents/Notre Dame/ND paper 2/Code/heatmap.png",base_height = 8,base_width = 8)
"""
