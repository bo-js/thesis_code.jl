using thesis_code, Random, Optimization, OptimizationNLopt, DelimitedFiles

T = 9000
burn = 1000

params = params_default()

grid = grids(; Nz = 11)
Z = grid[:Z]
X = grid[:X]
Y = grid[:Y]

Random.seed!(2077)

draw = rand(T+burn)

S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; 
    MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)

function HoardingCrit(thresholdsZ, params)
    thresholdszxy = [thresholdsZ[w] for w in 1:length(Z), i in 1:length(X), j in 1:length(Y)]
    crit = - optCrit(thresholdszxy, init; grid = grid, params = params, draw = draw, T = T, burn = burn, unconstr = true)
    return crit
end
x0 = zeros(Float64, length(Z))
func = Optimization.OptimizationFunction(HoardingCrit)

lower = zeros(Float64, length(Z))
lower .= -10.0
upper = zeros(Float64, length(Z))
upper .= 10.0

optHoarding = Optimization.OptimizationProblem(func, x0, params, lb = lower, ub = upper)

lossvals_hoarding = zeros(Float64, 0)

function callbackhoarding(state, loss_val)
    display(loss_val)
    push!(lossvals_hoarding, loss_val)
    return false
end

sol_hoarding = solve(optHoarding, NLopt.G_MLSL_LDS(); local_method = NLopt.LN_SBPLX(), callback = callbackhoarding, maxeval = 3000)
writedlm("output/opthoarding.txt", sol_hoarding.u)

