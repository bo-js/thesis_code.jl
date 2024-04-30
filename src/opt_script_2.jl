using thesis_code, Random, Optimization, OptimizationNLopt

T = 9000
burn = 1000

params = params_default()
grid = grids()

Random.seed!(2077)

draw = rand(T+burn)

x0 = [0.0, 0.0]

S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; 
    MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)

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

sol = solve(prob, NLopt.G_MLSL_LDS(); local_method = NLopt.LN_SBPLX(), callback = callback)


# model = Model(NLopt.Optimizer)
# set_optimizer_attribute(model, "algorithm", :LN_COBYLA)

# @variable(model, g1)
# @variable(model, g2)
# @variable(model, g3)
# @variable(model, g4)
# @variable(model, ui)
# @variable(model, 1 >= τ >= 0)


# @NLobjective(model, Max, optCrit(g1, g2, g3, g4, ui, τ, init; T = T, burn = burn, draw = draw))

# set_start_value(g1, 0)
# set_start_value(g2, 0)
# set_start_value(g3, 0)
# set_start_value(g4, 0)
# set_start_value(ui, 0)
# set_start_value(τ, 0)

# JuMP.optimize!(model)

# println("got ", objective_value(model), " at ", [value(g1), value(g2), value(g3), value(g4), value(ui), value(τ)])

# opt = Opt(:LN_COBYLA, 6)

# fcrit = pol ->  - optCrit(pol[1], pol[2], pol[3], pol[4], pol[5], pol[6], init; T = T, burn = burn, draw = draw, unconstr = false).c

# fconstr = pol -> optCrit(pol[1], pol[2], pol[3], pol[4], pol[5], pol[6], init; T = T, burn = burn, draw = draw, unconstr = false).d

# opt.lower_bounds = [-Inf, -Inf, -Inf, -Inf, -Inf, 0.]
# opt.upper_bounds = [Inf, Inf, Inf, Inf, Inf, 1]
# opt.xtol_rel = 1e-4

# opt.min_objective = fcrit
# inequality_constraint!(opt, fconstr, 1e-8)

# (minf, minx, ret) = NLopt.optimize(opt, x0)
# numevals = opt.numevals
# println("got $minf at $minx after $numevals iterations (returned $ret)")



