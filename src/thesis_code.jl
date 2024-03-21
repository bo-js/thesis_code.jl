module thesis_code

using LinearAlgebra

function valadd(Z::Vector, X::Vector, Y::Vector; p::Vector = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    return [z*(p ⋅ [1, x, y, x^2, y^2, x * y]) for z in Z, x in X, y in Y]
end

function homeprod(X::Vector, Y::Vector; p::Vector = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])    
    return [0.7 * maximum(valadd([x], Y, [1]; p = p)) for x in X]
end

function flow_surplus(Z, X, Y; p::Vector = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    va = valadd(Z, X, Y; p = p)
    b = homeprod(X, Y; p = p)
    return [va[w, i, j] - b[i] for  w in 1:length(Z), i in 1:length(X), j in 1:length(Y)]
end

include("surplus_vfi.jl")

export valadd
export homeprod
export flow_surplus
export SurplusVFI
export ES

include("helper_functions_dynamics.jl")


end