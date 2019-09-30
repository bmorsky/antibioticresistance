# Optimal control for the mean field model. Finds the optimal protocol for an
# initial condition of 900 susceptible bacteria by solving the state and adjoint
# equations. However, it does not prevent the emergence and establishment of the
# resistant type. Plots the time series of susceptible and resistance bacteria
# and the protocol.
using Plots, DiffEqBiological, Distributions, DifferentialEquations

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
μₐ = 0.000015
μ = 0.0000015
ζ = 1
η = 0.000001

# antibiotic off
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

p = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, μ, η)
prob = DiscreteProblem([500,500],(0.0,100.0),p)
jump_prob = JumpProblem(prob,Direct(),LV_antiOff)
sol = solve(jump_prob,FunctionMap())
plot(sol)

# antibiotic off
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
	μ₁, a + n₁ → a + n₂
	μ₂, n₂ → n₁
	#ζ, a → 0
	η, n₁ + n₂ → 2n₂
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ₁ μ₂ ζ η

p = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ₁, μ₂, ζ, η)
prob = DiscreteProblem([1000,0,10],(0.1,1000.0),p)
jump_prob = JumpProblem(prob,Direct(),LV_antiOn)
sol = solve(jump_prob,FunctionMap())
plot(sol)

# general reactions
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
	μ₁, a + n₁ → a + n₂
	μ₂, n₂ → n₁
	ζ, a → 0
	η, n₁ + n₂ → 2n₂
end b₁ d₁ b₂ d₂ γ₁₁ γ₁₂ γ₂₁ γ₂₂ α₁ α₂ μ₁ μ₂ ζ η

p = (b₁, d₁, b₂, d₂, γ₁₁, γ₁₂, γ₂₁, γ₂₂, α₁, α₂, μ₁, μ₂, ζ, η)
prob = DiscreteProblem([500,10,20],(0.0,100.0),p)
jump_prob = JumpProblem(prob,Direct(),LV_general)
sol = solve(jump_prob,FunctionMap())
plot(sol)
