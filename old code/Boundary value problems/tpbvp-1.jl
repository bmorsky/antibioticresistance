# Optimal control for the mean field model. Finds the optimal protocol for an
# initial condition of 900 susceptible bacteria by solving the state and adjoint
# equations. However, it does not prevent the emergence and establishment of the
# resistant type. Plots the time series of susceptible and resistance bacteria
# and the protocol.
using DifferentialEquations

#Parameters
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
init_pop = 5e3 #initial population size
final_time = 10000

function H(z)
	if z < 0
		return 1
	else
		return 0
	end
end

#System of ODEs: u[1] susceptible type, u[2] resistant type, u[3] adjoit equation for u[1], adjoint equation for u[2], control u[5].
function f(du,u,p,t)
	du[1] = (β₁-δ₁-μ)*u[1] - γ₁₁*u[1]^2 - (γ₁₂+η)*u[1]*u[2] + μ*u[2] - (α₁+μₐ)*u[1]*H(-u[3]*(α₁+μₐ)*u[1] + u[4]*(μₐ*u[1]-α₂*u[2]))
	du[2] = (β₂-δ₂-μ)*u[2] - γ₂₂*u[2]^2 + (η-γ₂₁)*u[1]*u[2] + μ*u[1] + (μₐ*u[1] - α₂*u[2])*H(-u[3]*(α₁+μₐ)*u[1] + u[4]*(μₐ*u[1]-α₂*u[2]))
	du[3] = -u[3]*(β₁ - δ₁ - μ - 2*γ₁₁*u[1] - (γ₁₂+η)*u[2] - (α₁+μₐ)*H(-u[3]*(α₁+μₐ)*u[1] + u[4]*(μₐ*u[1]-α₂*u[2]))) - u[4]*((η-γ₂₁)*u[2] + μ + μₐ*H(-u[3]*(α₁+μₐ)*u[1] + u[4]*(μₐ*u[1]-α₂*u[2])))
	du[4] = -u[3]*(-(γ₁₂+η)*u[1] + μ) - u[4]*(β₂ - δ₂ - μ - 2*γ₂₂*u[2] + (η-γ₂₁)*u[1] - α₂*H(-u[3]*(α₁+μₐ)*u[1] - u[4]*(μₐ*u[1]-α₂*u[2])))
end

#Solve the system
function bc!(residual, u, p, t)
	residual[1] = u[1][1] - init_pop # the solution at the middle of the time span should be -pi/2
	residual[2] = u[1][2] # the solution at the end of the time span should be pi/2
    residual[3] = u[end][3] - 1 # the solution at the middle of the time span should be -pi/2
    residual[4] = u[end][4] - 1 # the solution at the end of the time span should be pi/2
end
u0 = [init_pop,0,1,1] #initial conditions: init_pop susceptible, 0 resistance, adjoint equations 1 and 1 for control minimizing the total number of bacteria at the final time, and the antibiotic on i.e. a = 1
tspan = (0.0,final_time) #time span
prob = TwoPointBVProblem(f,bc!,u0,tspan) #the problem to solve
sol = solve(prob,MIRK4(),dt = 0.5) #solve
#sol = solve(prob,Tsit5()) #solve
#plot(sol) #plot in julia

#for plotting in R
output = [sol[1,:] sol.t[:] zeros(length(sol[1,:]),1); sol[2,:] sol.t[:] ones(length(sol[2,:]),1)] #output
protocol = [sign.(-(α₁+μₐ)*sol[3,:].*sol[1,:]  + sol[4,:].*(μₐ*sol[1,:]-α₂*sol[2,:])) sol.t[:]] #the protocol
#R code to generate figure
using RCall
@rput output protocol;
R"""
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

output <- as.data.frame(output)
names(output) <- c("Number","Time","Type")
protocol <- as.data.frame(protocol)
names(protocol) <- c("Antibiotic","Time")

ts <- ggplot() + scale_x_continuous(expand = c(0, 0))

output <- ts + geom_line(data=output,aes(x=Time,y=Number,group=Type),size=1) + ggtitle("Time series for the optimal control")+ scale_color_manual(values=c("blue","red"),labels = c("Susceptible", "Resistant")) + aes(color = factor(Type)) + theme(legend.title=element_blank(),legend.position=c(.65,.85)) + labs(x = "Time", y = "Bacterial load") + scale_y_continuous(expand = c(0, 0))
protocol <- ts + geom_line(data=protocol,aes(x=Time,y=Antibiotic),size=1) + ggtitle("Time series for the optimal control") + theme(legend.title=element_blank()) + labs(x = "Time", y = "Antibiotic dose") + scale_y_continuous(expand = c(0, 0))

plot_out <- plot_grid(output,protocol,labels=letters[1:2],ncol=2)
save_plot(plot_out,filename="/home/bmorsky/antibiotic/tpbvp.png",base_height = 5,base_width = 10)
"""
