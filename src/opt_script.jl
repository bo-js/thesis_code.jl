using thesis_code, Optim

T = 10000
burn = 2000

params = params_default()
grid = grids()

draw = rand(T+burn)

x0 = [0.0, 0.0]
f = pol -> -optCrit(pol[1], pol[2]; T = T, burn = burn, draw = draw, unconstr = true)
res = optimize(f, x0, NelderMead())

Optim.minimizer(res)