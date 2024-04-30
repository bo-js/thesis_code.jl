using Distributed

addprocs(2)

@everywhere using thesis_code, BlackBoxOptim, Random

@everywhere T = 20000
@everywhere burn = 2000

@everywhere params = params_default()
@everywhere grid = grids()

@everywhere Random.seed!(2077)

@everywhere draw = rand(T+burn) 

@everywhere S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
@everywhere init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; 
    MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)

#@everywhere origfitness(x::Vector{Float64}) =  - (optCrit(0,0,0,0,x[1], x[2], init;  params = params, T = T, burn = burn, draw = draw, unconstr = false).c)

@everywhere Dim = 2
@everywhere lowbounds = [-2.0, 0.0]
@everywhere highbounds = [2.0, 1.0]

@everywhere SearchRange = [(lowbounds[i], highbounds[i]) for i in 1:Dim]

# @everywhere constraint(x::Vector{Float64}) = max(0, optCrit(0,0,0,0,x[1], x[2], init; params = params, T = T, burn = burn, draw = draw, unconstr = false).d)

@everywhere K = 1e3

@everywhere function penalized_fitness(x)
    out = optCrit(0,0,0,0, x[1], x[2], init; params = params, T=T, burn = burn, draw = draw, unconstr = false)
    return - out.c + K * max(0, out.d)
end

opt = bbsetup(penalized_fitness; SearchRange = SearchRange, Workers = workers(), MaxFuncEvals = 5000, Method=:adaptive_de_rand_1_bin_radiuslimited)

res = bboptimize(opt, [0.0, 0.0])

best = best_candidate(res)