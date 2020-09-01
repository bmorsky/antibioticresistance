# Simulations for the stochastic model for a variety of protocols. Results plotted as heatmaps.
using Plots, Distributions

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

maxOff = 20; #maximum duration the antibiotic is off
maxOn = 20; #maximum duration the antibiotic is on
avg_work = zeros(maxOn*(maxOff+1),3);
std_work = zeros(maxOn*(maxOff+1),3);
avg_t = zeros(maxOn*(maxOff+1),3);
std_t = zeros(maxOn*(maxOff+1),3);

#Dynamics for when the antibiotic is on
function antibioticOn(period,output)
	N = zeros(period+1,2);
	N[1,1] = output[1];
	N[1,2] = output[2];
	for t = 2:period+1
		N[t,1] = rand(Poisson(max(d*N[t-1,1]*(1 - dot(N[t,:],Ia)/K),0)));
		N[t,2] = rand(Poisson(max(w_r*N[t-1,2]*(1 - dot(N[t,:],Ia_r)/K),0)));
		if N[t,2] > 1000
			N[t,2] = 1000;
		end
		if N[t,1] > 0
			mutants = rand(Binomial(N[t,1], m));
			N[t,1] = N[t,1] - mutants;
			N[t,2] = N[t,2] + mutants;
		end
	end
	return N[2:end,:]
end

#Dynamics for when the antibiotic is off
function antibioticOff(period,output)
	N = zeros(period+1,2);
	N[1,1] = output[1];
	N[1,2] = output[2];
	for t = 2:period+1
		N[t,1] = rand(Poisson(max((w+s)*N[t-1,1]*(1 - dot(N[t,:],I)/K),0)));
		N[t,2] = rand(Poisson(max(w*N[t-1,2]*(1 - dot(N[t,:],I_r)/K),0)));
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

#Runs simulations
n = zeros(1,2);
count = 1;
for dur_on = 5:5:5*maxOn
	for dur_off = 0:5:5*maxOff
		if dur_on == 0 && dur_off == 0
			break;
		end
		work = zeros(num_sims,2);
		for sim = 1:num_sims
			#Initialize population
			n = zeros(1,2);
			n[1,1] = init_pop; #initial number of susceptible type for the protocol
			n[1,2] = 0; #initial number of resistant type for the protocol
			t = 0;
			while sum(n[end,:]) > 0 && t <= 100 #n[end,2] < 2000 && n[end,1] < 2000
				n = vcat(n, antibioticOn(dur_on,n[end,:]));
				n = vcat(n, antibioticOff(dur_off,n[end,:]));
				t += dur_on + dur_off;
			end
			work[sim,:] = [convert(Float64, sum(n[end,:]) == 0.0) t];
		end
		avg_work[count,:] = [dur_on dur_off mean(work[:,1])]
		std_work[count,:] = [dur_on dur_off std(work[:,1])]
		avg_t[count,:] = [dur_on dur_off mean(work[:,2])]
		std_t[count,:] = [dur_on dur_off std(work[:,2])]
		count += 1;
	end
end

#R code to generate heatmaps for the probability of success of the protocols, the average time to success, and the variances of these
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
