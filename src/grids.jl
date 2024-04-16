
function grids(; Nx = 21, Ny = 21, Nz = 51, β1 = 2.148, β2 = 12.001, ρ = 0.99, σ = 0.071)

     trim = 0.01

     X = collect(LinRange(trim, 1 - trim, Nx))
     Y = collect(LinRange(trim, 1 - trim, Ny))
     a = collect(LinRange(trim, 1-trim, Nz))

     l = pdf(Beta(β1, β2), X)
     l = l ./ sum(l)

     Σ = [1 ρ
          ρ 1]

     Z = exp.(σ * quantile(Normal(), a))

     cop = GaussianCopula(Σ)

     P = [pdf(cop, [a[i], a[j]]) for i in 1:Nz, j in 1:Nz]
     Π = [P[i, j]/sum(P[i, :]) for i in 1:Nz, j in 1:Nz]
     
     return Dict(:X => X, :Y => Y, :Z => Z, :Π => Π, :l => l)
end

export grids

