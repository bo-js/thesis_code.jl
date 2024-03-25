function u_plus(uxt::Matrix, hxt::Array, S::Array, statet::Vector, t::Number; δ::Number = δ)
    return [uxt[t, i] + sum(hxt[t, i, j] * ((S[statet[t], i, j] < 0) + δ * (S[statet[t], i, j] ≥ 0)) for j in 1:length(S[1, 1, :])) for i in 1:length(S[1, :, 1])] 
end

export u_plus

function h_plus(hxt::Array, S::Array, statet::Vector, t::Number; δ::Number = δ)
    return [(1 - δ) * (S[statet[t], i, j]) * hxt[t, i, j] for i in 1:length(S[1, :, 1]), j in 1:length(S[1, 1, :])]
end

export h_plus

function LT(uxt::Matrix, hxt::Array, t::Number; s::Number = s)
    u_plus = u_plus(uxt, hxt, S, statet, t)
    h_plus = h_plus(hxt, S, statet, t)
    return sum(u_plus[i] for i in 1:length(uxt[t, :])) + s * sum(h_plus[i, j] for i in 1:length(h[t, :, 1]) j in 1:length(h[t, 1, :]))
end

export LT

function J(uxt::Matrix, hxt::Array, S::Array, statet::Vector, t::Number; s::Number = s)
    u_plus = u_plus(uxt, hxt, S, statet, t)
    h_plus = h_plus(hxt, S, statet, t)
    L = LT(uxt, hxt, t)

    return [
        sum(u_plus[i] * max(S[statet[t], i, j], 0)/L for i in 1:length(S[1, :, 1])) + 
        sum(s * h_plus[i, jprime] * max(S[statet[t], i, j] - S[statet[t], i, jprime], 0)/L for i in 1:length(S[1, :, 1]) jprime in 1:length(S[1, 1, :]))
    ]
    
end

export J

function tightness(uxt::Matrix, hxt::Array, S::Array, statet::Vector, t::Number; α::Number = α, c0::Number = c0, c1::Number = c1, ω::Number = ω)
    J = J(uxt, hxt, S, statet, t)
    return (sum((α * J[j]/c0)^(1/c1) for j in 1:Ny)/LT(uxt, hxt, t))^(c1/(c1 + ω))
end

export tightness

function vacy(uxt::Matrix, hxt::Array, S::Array, statet::Vector, t::Number; α::Number = α, c0::Number = c0, c1 = c1, ω = ω)
    
    return ((α/(c0 * tightness(uxt, hxt, S, statet, t)^ω)) .* J(uxt, hxt, S, statet, t)).^(1/c1)

end

export vacy

function uh_next(uxt::Matrix, hxt::Array, S::Array, statet::Vector, t::Number; α::Number = α, ω::Number = ω, s::Number = s)

    u_plus = u_plus(uxt, hxt, S, statet, t)
    h_plus = u_plus(hxt, S, statet, t)
    L = LT(uxt, hxt, t)
    v_y = vacy(uxt, hxt, S, statet, t)
    V = sum(v_y)
    λ = min(α * L^ω * V^(1 - ω), L, V)/L

    u_next = [u_plus[i] * (1 - sum(λ * (v[j]/V)*(S[statet[t], i, j] > 0) for j in 1:Ny)) for i in 1:Nx]


    h_next = [
        h_plus[i, j] * 
        
        (1 - sum(s * λ * (v_y[j_prime]/V) * (S[statet[t], i, j_prime] > S[statet[t], i, j]) for j_prime in 1:Ny)
        
        + sum(h_plus[i, j_prime] * s * λ * (v_y[j]/V) * (S[statet[t], i, j] > S[statet[t], i, j_prime]) for j_prime in 1:Ny)
        
        + u_plus[i] * λ * (v_y[j]/V) * (S[statet[t], i, j] ≥ 0)) 
        
        for i in 1:Nx, j in 1:Ny
    ]

    next::NamedTuple = (u = u_next, h = h_next)
    
    return next

end

export uh_next

# NEED TO FIGURE OUT HOW TO GET STATET IN - POTENTIALLY DONE IN LOOP
# function uh_dynamics_sim(u0::Vector, h0::Matrix, S::Array, statet::Vector, T::Integer)
#     uxt = zeros(T, Nx)
#     hxt = zeros(T, Nx, Ny)

#     uxt[1, :] = u0
#     hxt[1, :, :] = h0

#     for t in 2:T
#         next = uh_next(uxt, hxt, S, statet, t - 1)
#         uxt[t, :] = next.u
#         hxt[t, :, :] = next.h
#     end

#     return Dict(uxt => :uxt, hxt => :hxt)
# end

# export uh_dynamics_sim