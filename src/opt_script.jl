using thesis_code, Optim

T = 50000
burn = 5000

params = params_default()
grid = grids()

draw = rand(T+burn)

x0 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.2]
f = pol -> -optCrit(pol[1], pol[2], pol[3], pol[4], pol[5], pol[6]; T = T, burn = burn, draw = draw, unconstr = true)
res = optimize(f, x0, NelderMead())

argmin = Optim.minimizer(res)

