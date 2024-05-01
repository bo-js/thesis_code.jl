## SETUP
using thesis_code, Random, Optimization, OptimizationNLopt, DelimitedFiles

T = 9000
burn = 1000

params = params_default()
grid = grids()

Random.seed!(2077)

draw = rand(T+burn)

S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; 
    MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)


## UI only
x0 = [0.0,0.0]
f = OptimizationFunction((pol, params) -> - optCrit(pol[1], pol[2], init; params = params, T = T, burn = burn, draw = draw, unconstr = true))
prob = Optimization.OptimizationProblem(f, x0, params, lb = [-1.0, 0.0], ub = [1.0, 1.0])

lossvals = zeros(0)
ui_guesses = zeros(0)
tax_guesses = zeros(0)

function callback(state, loss_val)
    display(loss_val)
    push!(ui_guesses, state[1])
    push!(tax_guesses, state[2])
    push!(lossvals, loss_val)
    return false
end

sol = solve(prob, NLopt.G_MLSL_LDS(); local_method = NLopt.LN_SBPLX(), callback = callback, maxeval = 3000)
writedlm("output/optui.txt",sol.u)
## UI with Threshold
x0 = [0.0,0.0, 0.0]

f2 = OptimizationFunction((pol, params) -> - optCrit(pol[1], pol[2], init;threshold = pol[3], params = params, T = T, burn = burn, draw = draw, unconstr = true))
prob2 = Optimization.OptimizationProblem(f2, x0, params, lb = [-1.0, 0.0, -5.0], ub = [1.0, 1.0, 1.0])

lossvals2 = zeros(0)
ui_guesses2 = zeros(0)
tax_guesses2 = zeros(0)
threshold_guesses = zeros(0)

function callback2(state, loss_val)
    display(loss_val)
    push!(ui_guesses, state[1])
    push!(tax_guesses, state[2])
    push!(threshold_guesses, state[3])
    push!(lossvals2, loss_val)
    return false
end

sol2 = solve(prob2, NLopt.G_MLSL_LDS(); local_method = NLopt.LN_SBPLX(), callback = callback2, maxeval = 3000)
writedlm("output/optui_threshold.txt", sol2.u)

## Subsidy Only
x0 = [0.0,0.0,0.0,0.0,0.0]

f3 = OptimizationFunction((pol, params) -> - optCrit(pol[1], pol[2], pol[3], pol[4], 0, pol[5], init; params = params, T = T, burn = burn, draw = draw, unconstr = true))
prob3 = Optimization.OptimizationProblem(f3, x0, params, lb = [-10.0, -5.0, -5.0, -5.0, 0.0], ub = [10.0, 5.0, 5.0, 5.0, 1.0])

lossvals3 = zeros(0)
pol1_guesses = zeros(0)
pol2_guesses = zeros(0)
pol3_guesses = zeros(0)
pol4_guesses = zeros(0)
tax_guesses3 = zeros(0)

function callback3(state, loss_val)
    display(loss_val)
    push!(pol1_guesses, state[1])
    push!(pol2_guesses, state[2])
    push!(pol3_guesses, state[3])
    push!(pol4_guesses, state[4])
    push!(tax_guesses3, state[5])
    push!(lossvals3, loss_val)
    return false
end

sol3 = solve(prob3, NLopt.G_MLSL_LDS(); local_method = NLopt.LN_SBPLX(), callback = callback3, maxeval = 3000)
writedlm("output/optsub.txt", sol3)

## Subsity and UI
x0 = [0.0,0.0,0.0,0.0,0.0,0.0]

f4 = OptimizationFunction((pol, params) -> - optCrit(pol[1], pol[2], pol[3], pol[4], pol[5], pol[6], init; params = params, T = T, burn = burn, draw = draw, unconstr = true))
prob4 = Optimization.OptimizationProblem(f4, x0, params, lb = [-10.0, -5.0, -5.0, -5.0, -1.0, 0.0], ub = [10.0, 5.0, 5.0, 5.0, 1.0, 1.0])

lossvals4 = zeros(0)
pol1_guesses_full = zeros(0)
pol2_guesses_full = zeros(0)
pol3_guesses_full = zeros(0)
pol4_guesses_full = zeros(0)
ui_guesses4 = zeros(0)
tax_guesses4 = zeros(0)

function callback4(state, loss_val)
    display(loss_val)
    push!(pol1_guesses_full, state[1])
    push!(pol2_guesses_full, state[2])
    push!(pol3_guesses_full, state[3])
    push!(pol4_guesses_full, state[4])
    push!(ui_guesses4, state[5])
    push!(tax_guesses4, state[6])
    push!(lossvals4, loss_val)
    return false
end

sol4 = solve(prob4, NLopt.G_MLSL_LDS(); local_method = NLopt.LN_SBPLX(), callback = callback4, maxeval = 3000)
writedlm("output/optsub_ui.txt", sol4.u)