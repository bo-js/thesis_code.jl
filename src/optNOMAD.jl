using thesis_code, Random, NOMAD

T = 10000
burn = 1000

params = params_default()
grid = grids()

Random.seed!(2077)

draw = rand(T+burn)

S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; 
    MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)

function eval_fct(pol)
    sim = optCrit(0,0,0,0, pol[1], pol[2], init; params = params, T = T, burn = burn, draw = draw)
    bb_outputs = [sim.c, sim.d]
    success = (bb_outputs[1] != NaN)
    count_eval = true
    return (success, count_eval, bb_outputs)
end

pb = NomadProblem(
    2,
    2,
    ["OBJ", "EB"],
    eval_fct;
    lower_bound = [-2.0, 0.0],
    upper_bound = [2.0, 1.0]
)

result = solve(pb, [0.0, 0.0])