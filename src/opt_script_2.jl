using Optimization, OptimizationNLopt, thesis_code, Random

T = 50000
burn = 5000

params = params_default()
grid = grids()

Random.seed!(2077)

draw = rand(T+burn)

x0 = [0.0, 0.0]

f = OptimizationFunction((pol, params) -> - optCrit(0, 0, 0, 0, pol[1], pol[2]; params = params, T = T, burn = burn, draw = draw, unconstr = true))
prob = Optimization.OptimizationProblem(f, x0, params, lb = [-Inf, 0.], ub = [ Inf, 1])

sol = solve(prob, NLopt.LN_COBYLA())


# model = Model(NLopt.Optimizer)
# set_optimizer_attribute(model, "algorithm", :LN_COBYLA)

# @variable(model, g1)
# @variable(model, g2)
# @variable(model, g3)
# @variable(model, g4)
# @variable(model, ui)
# @variable(model, 1 >= τ >= 0)

# fcrit = (x1, x2, x3, x4, ui, τ) -> optCrit(x1, x2, x3, x4, ui, τ; T = T, burn = burn, draw = draw, unconstr = false).c

# register(model, :fcrit, 6, fcrit, autodiff=true)

# @NLobjective(model, Max, optCrit(g1, g2, g3, g4, ui, τ; T = T, burn = burn, draw = draw))

# set_start_value(g1, 0)
# set_start_value(g2, 0)
# set_start_value(g3, 0)
# set_start_value(g4, 0)
# set_start_value(ui, 0)
# set_start_value(τ, 0)

# JuMP.optimize!(model)

# println("got ", objective_value(model), " at ", [value(g1), value(g2), value(g3), value(g4), value(ui), value(τ)])

# opt = Opt(:LN_COBYLA, 6)

# function fcrit(pol::Vector) 
#     crit =  - optCrit(pol[1], pol[2], pol[3], pol[4], pol[5], pol[6]; T = T, burn = burn, draw = draw, unconstr = false).c
#     return crit
# end

# function fconstr(pol::Vector) 
#     def = optCrit(pol[1], pol[2], pol[3], pol[4], pol[5], pol[6]; T = T, burn = burn, draw = draw, unconstr = false).d
#     return def
# end

# opt.lower_bounds = [-Inf, -Inf, -Inf, -Inf, -Inf, 0.]
# opt.upper_bounds = [Inf, Inf, Inf, Inf, Inf, 1]
# opt.xtol_rel = 1e-4

# opt.min_objective = fcrit
# inequality_constraint!(opt, fconstr, 1e-8)

# (minf, minx, ret) = NLopt.optimize(opt, x0)
# numevals = opt.numevals
# println("got $minf at $minx after $numevals iterations (returned $ret)")



