using Distributions
using BivariateCopulas

Nx = 21
Ny = 21
Nz = 51

trim = 0.01

X = collect(LinRange(trim, 1 - trim, Nx))
Y = collect(LinRange(trim, 1 - trim, Ny))
a = collect(LinRange(trim, 1-trim, Nz))

σ = 0.071
ρ = 0.92 # FOR SOME REASON NOT CURRENTLY WORKING WITH PAPER VALUE OF 0.999

Z = exp.(σ * quantile(Normal(), a))

cop = Gaussian(ρ)

P = [density(cop, a[i], a[j]) for i in 1:Nz, j in 1:Nz]
Π = [P[i, j]/sum(P[i, :]) for i in 1:Nz, j in 1:Nz]
