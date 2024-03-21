function u_plus(uxt::Matrix, hxt::Array, S::Array, statet::Vector, t::Number; δ::Number = 0.013)
    return [uxt[t, i] + sum(hxt[t, i, j] * ((S[statet[t], i, j] < 0) + δ * (S[statet[t], i, j] ≥ 0)) for j in 1:length(S[1, 1, :])) for i in 1:length(S[1, :, 1])] 
end

export u_plus

function h_plus(hxt::Array, S::Array, statet::Vector, t::Number; δ::Number = 0.013)
    return [(1 - δ) * (S[statet[t], i, j]) * hxt[t, i, j] for i in 1:length(S[1, :, 1]), j in 1:length(S[1, 1, :])]
end

export h_plus

function LT(uxt::Matrix, hxt::Array, t::Number; s::Number = 0.027)
    return sum(uxt[t, i] for i in 1:length(uxt[t, :])) + s * sum(hxt[t, i, j] for i in 1:length(h[t, :, 1]) j in 1:length(h[t, 1, :]))
end

export LT

function J(uxt, hxt, S::Array, statet::Vector, t::Number; s::Number = 0.027, δ::Number = 0.013)
    u_plus = u_plus(uxt, hxt, S, statet, t; δ = δ)
    h_plus = h_plus(hxt, S, statet, t; δ)
    L = LT(uxt, hxt, t; s = s)

    return [
        sum(u_plus[i] * max(S[statet[t], i, j], 0)/L for i in 1:length(S[1, :, 1])) + 
        sum(s * h_plus[i, jprime] * max(S[statet[t], i, j] - S[statet[t], i, jprime], 0)/L for i in 1:length(S[1, :, 1]) jprime in 1:length(S[1, 1, :]))
    ]
    
end

export J

function tightness(uxt, hxt, S, statet, t; s::Number = 0.027, α::Number = 0.497, δ::Number =  0.013, c0::Number = 0.028, c1 = 0.084, ω = 1/2)
    J = J(uxt, hxt, S, statet, t; s = s, δ = δ)
    return (sum((α * J[j]/c0)^(1/c1) for j in 1:Ny)/LT(uxt, hxt, t; s = s))^(c1/(c1 + ω))
end

export tightness

function vacy(uxt, hxt, S, statet, t; s::Number = 0.027, α::Number = 0.497, δ::Number =  0.013, c0::Number = 0.028, c1 = 0.084, ω = 1/2)
    
    return ((α/(c0 * tightness(uxt, hxt, S, statet, t; α = α, δ = δ, c0 = c0, c1 = c1, ω = ω)^ω)) .* J(uxt, hxt, S, statet, t; s = s, δ = δ)).^(1/c1)

end

export vacy
