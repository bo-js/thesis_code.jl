using thesis_code
using Distributions
using Copulas

# Parameters
## Interest Rate
r = 0.05/52
## Matching
α = 0.497
ω = 1/2
## Search Intensity
s = 0.027
## Vacancy Costs
c0 = 0.028
c1 = 0.084
## Exogeneous Separation
δ = 0.013
## Productivity Shocks
σ = 0.071
ρ = 0.999
## Worker Heterogeneity
β1 = 2.148
β2 = 12.001
## Value Added
p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596]

# Define Grids, Distributions of Workers, and Productivity Markov Process
grd = grids(; β1 = β1, β2 = β2, ρ = ρ, σ = σ)
X = grd[:X]
Y = grd[:Y]
Z = grd[:Z]
Π = grd[:Π]
l = grd[:l]

# Find Surplus Function
S = SurplusVFI(Z, X, Y, Π; δ = δ, r = r)

# Load Data
prod = parse.(Float64, readlines("/Users/bojs/Desktop/LR_changes/1_data/data_for_fortran/output/gdp_w.txt"))

# Define Initial Conditions
init = H0iter(Z, X, Y, S, l; δ = δ, s=s, α = α, c0 = c0, c1 = c1, ω = ω)


# Simulate
sim = uh_dynamics_sim(prod, init.u, init.h, S, l; Z = Z, X = X, Y = Y, δ = δ, s = s, α = α, c0 = c0, c1 = c1, ω = ω, p = p)
uxt = sim[:uxt]
hxyt = sim[:hxyt]
statet = sim[:statet]
GDP_zt = sim[:GDP_zt]
