using thesis_code, PrettyTables, Latexify, Random, LinearAlgebra


T = 9000
burn = 1000

Random.seed!(2077)
draw = rand(T+burn)

params = params_default()
grid = grids()

Random.seed!(2077)

draw = rand(T+burn)

# Get Params
p = params.p
c0 = params.c0
c1 = params.c1

# Define Grids 
X = grid[:X]
Y = grid[:Y]
Z = grid[:Z]

# Steady State
S0 = SurplusVFI(grid[:Z], grid[:X], grid[:Y], grid[:Π]; δ = params.δ, r = 0.05/52, p = params.p)
init = H0iter(grid[:Z], grid[:X], grid[:Y], S0, grid[:l]; MaxIter = 5000, tol = 10e-5, δ = params.δ, s = params.s, α = params.α, c0 = params.c0, c1 = params.c1, ω = params.ω)
h0 = init.h
u0 = init.u

# Value added and home production
va = valadd(Z, X, Y; p = p)
b = homeprod(X, Y; p = p)


# Baseline
fs_bsline = flow_surplus(va, b)

sim_bsline = sim_draws(grid, params, fs_bsline, init; T = T, burn = burn, draw = draw)

swf_emp_bsline = sum(sum(sim_bsline.hxyt[t, :, :] .* (va[sim_bsline.statet[t], :, :])) for t in (burn+1):(T+burn))/T
swf_unemp_bsline = sum(dot(sim_bsline.uxt[t, :], b) for t in (burn+1):(T+burn))/T
swf_vaccost_bsline = sum(- (c0/(1 + c1)) * sum(sim_bsline.vacyt[t, :] .^ (1 + c1)) for t in (burn+1):(T+burn))/T

# Optimal Subsidy
g = parse.(Float64, readlines("output/optsub.txt"))

sub = [(Z[w] < g[1]) * g[2] for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]

tax = [ - (g[3] * va[w, i, j]) for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]

fs_sub = flow_surplus(va + sub + tax, b)

sim_sub = sim_draws(grid, params, fs_sub, init; T = T, burn = burn, draw = draw)

swf_emp_sub = sum(sum(sim_sub.hxyt[t, :, :] .* va[sim_sub.statet[t], :, :]) for t in (burn+1):(T+burn))/T
swf_subsidy = sum(sum(sim_sub.hxyt[t, :, :] .* sub[sim_sub.statet[t], :, :]) for t in (burn+1):(T+burn))/T
swf_tax_sub = sum(sum(sim_sub.hxyt[t, :, :] .* tax[sim_sub.statet[t], :, :]) for t in (burn+1):(T+burn))/T
swf_unemp_sub = sum(dot(sim_sub.uxt[t, :], b) for t in (burn+1):(T+burn))/T
swf_vaccost_sub = sum(- (c0/(1 + c1)) * sum(sim_sub.vacyt[t, :] .^ (1 + c1)) for t in (burn+1):(T+burn))/T

# Optimal UI
g = parse.(Float64, readlines("output/optui.txt"))

uix = g[1] .* [sum(va[Int((length(Z)+1)/2), i, :] .* h0[i, :])/sum(h0[i, :]) for i in 1:length(X)]
tax_ui = [ - (g[2] * va[w, i, j]) for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]

fs_ui = flow_surplus((1 - g[2]) .* va, b.+ uix)

sim_ui = sim_draws(grid, params, fs_ui, init; T = T, burn = burn, draw = draw)

swf_emp_ui = sum(sum(sim_ui.hxyt[t, :, :] .* va[sim_ui.statet[t], :, :]) for t in (burn+1):(T+burn))/T
swf_tax_ui = sum(sum(sim_ui.hxyt[t, :, :] .* tax_ui[sim_ui.statet[t], :, :]) for t in (burn+1):(T+burn))/T
swf_unemp_ui = sum(dot(sim_ui.uxt[t, :], b) for t in (burn+1):(T+burn))/T
swf_uipayments = sum(dot(sim_ui.uxt[t, :], uix) for t in (burn+1):(T+burn))/T
swf_vaccost_ui = sum(- (c0/(1 + c1)) * sum(sim_ui.vacyt[t, :] .^ (1 + c1)) for t in (burn+1):(T+burn))/T

# # Both
# g = parse.(Float64, readlines("output/optsub_ui.txt"))

# sub = [(Z[w] < g[1]) * g[2] - (g[4] * va[w, i, j]) for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]

# uix = g[3] .* [sum(va[Int((length(Z)+1)/2), i, :] .* h0[i, :])/sum(h0[i, :]) for i in 1:length(X)]

# fs_both = flow_surplus(va .+ sub, b .+ uix)

# sim_both = sim_draws(grid, params, fs_both, init; T=T, burn = burn, draw = draw)

# swf_emp_both = sum(sum(sim_both.hxyt[t, :, :] .* (va[sim_both.statet[t], :, :] + sub[sim_both.statet[t], :, :])) for t in (burn+1):(T+burn))/T
# swf_unemp_both = sum(dot(sim_both.uxt[t, :], b .+ uix) for t in (burn+1):(T+burn))/T
# swf_vaccost_both = sum(- (c0/(1 + c1)) * sum(sim_both.vacyt[t, :] .^ (1 + c1)) for t in (burn+1):(T+burn))/T


# Table
header = LatexCell.(["Policy", latexify("Δ") * " Value Added", latexify("Δ") * " Home Prod.", latexify("Δ") * " Vacancy Costs", "UI Payments", "Subsidy Payments", "Tax", latexify("Δ") * " Social Welfare"])
cpol = ["UI", "Subsidy"]
cva = [swf_emp_ui - swf_emp_bsline, swf_emp_sub - swf_emp_bsline]
sub = [0, swf_subsidy]
tax = [swf_tax_ui, swf_tax_sub]
chome = [swf_unemp_ui - swf_unemp_bsline, swf_unemp_sub - swf_unemp_bsline]
ui = [swf_uipayments, 0]
cvac = [swf_vaccost_ui - swf_vaccost_bsline, swf_vaccost_sub - swf_vaccost_bsline]
cswf = [
    swf_emp_ui + swf_unemp_ui + swf_vaccost_ui - (swf_emp_bsline + swf_unemp_bsline + swf_vaccost_bsline),
    swf_emp_sub + swf_unemp_sub + swf_vaccost_sub - (swf_emp_bsline + swf_unemp_bsline + swf_vaccost_bsline),
]

tab = pretty_table(String,
    hcat(cpol, cva, chome, cvac, ui, sub, tax, cswf);
    header = header,
    backend = Val(:latex)
)

open("output/swf_table.tex", "w") do io
    print(io, tab)
end

