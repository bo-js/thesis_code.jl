using thesis_code

# Parameters
## Interest Rate
r = 0.05
## Matching
α = 0.097
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
ρ = 0.99
## Worker Heterogeneity
β1 = 2.148
β2 = 12.001
## Value Added
p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596]

# Define Grids, Distributions of Workers, and Productivity Markov Process
include("grids.jl")

# Find Surplus Function
S = SurplusVFI(Z, X, Y, Π; δ = δ, r = r)

# Load Data
prod = parse.(Float64, readlines("/Users/bojs/Desktop/LR_changes/1_data/data_for_fortran/output/gdp_w.txt"))
