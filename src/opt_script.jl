using thesis_code, Optim, Random

T = 9000
burn = 1000

params = params_default()
grid = grids()

Random.seed!(2077)
draw = rand(T+burn)

S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; 
    MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)


x0 = [0.24, 2.21e-5, 0.0]
f = pol -> -optCrit(pol[1],pol[2], init; threshold = pol[3], T = T, burn = burn, draw = draw, unconstr = true)
res = Optim.optimize(f, x0, NelderMead())

argmin = Optim.minimizer(res)

