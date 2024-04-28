using Distributed

addprocs(8)

@everywhere using thesis_code, BlackBoxOptim, Random

@everywhere T = 50000
@everywhere burn = 5000

@everywhere params = params_default()
@everywhere grid = grids()

@everywhere Random.seed!(2077)

@everywhere draw = rand(T+burn)

@everywhere origfitness(x::Vector{Float64}) =  - (optCrit(0,0,0,0,x[1], x[2];  params = params, T = T, burn = burn, draw = draw, unconstr = false).c)

@everywhere Dim = 2
@everywhere lowbounds = [-2.0, 0.0]
@everywhere highbounds = [2.0, 1.0]

@everywhere SearchRange = [(lowbounds[i], highbounds[i]) for i in 1:Dim]

@everywhere constraint(x::Vector{Float64}) = max(0, optCrit(0,0,0,0,x[1], x[2]; params = params, T = T, burn = burn, draw = draw, unconstr = false).d)

@everywhere K = 1e3

@everywhere penalized_fitness(x) = origfitness(x) + K * constraint(x)

opt = bbsetup(penalized_fitness; SearchRange = SearchRange, Workers = workers())

res = bboptimize(opt)