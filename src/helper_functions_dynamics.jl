function u_plus(ux::Vector, hxy::Matrix, Sxy::Matrix; δ::Number = δ)
   return ux .+ vec(sum(((Sxy .< 0) + δ.*(Sxy .≥ 0)) .* hxy , dims = 2))
end

export u_plus

function h_plus(hxy::Matrix, Sxy:: Matrix; δ::Number = δ)
    return (1 - δ) .* (Sxy .≥ 0) .* hxy
end

export h_plus

function LT(u_plus::Vector, h_plus::Matrix; s::Number = s)
    return sum(u_plus) + s * sum(h_plus)
end

export LT

function J(u_plus::Vector, h_plus::Matrix, L::Number, Sxy::Matrix; s::Number = s)
    
    return (sum(u_plus[i] .* max.(Sxy[i, :], 0)./L for i in 1:length(u_plus)) 
        
            .+ sum(s .* h_plus[i, jprime] .* max.(Sxy[i, :] .- Sxy[i, jprime], 0)./L 
                
                    for i in 1:length(u_plus), jprime in 1:length(h_plus[1, :])))
    
end

export J

function tightness(Jy::Vector, L::Number; α::Number = α, c0::Number = c0, c1::Number = c1, ω::Number = ω)
    
    return (sum((α * Jy[j]/c0)^(1/c1) for j in 1:length(Jy))/L)^(c1/(c1 + ω))

end

export tightness

function vacy(Jy::Vector, θ::Number; α::Number = α, c0::Number = c0, c1 = c1, ω = ω)
    
    return ((α/(c0 * θ^ω)) .* Jy).^(1/c1)

end

export vacy

function uh_next(u_plus::Vector, h_plus::Matrix, v_y::Vector, L::Number, Sxy::Matrix, l::Vector; α::Number = α, ω::Number = ω, s::Number = s)

    V = sum(v_y)
    λ = min(α * L^ω * V^(1 - ω), L, V)/L
    λvV = λ .* v_y ./V


    h_next = (
        h_plus .* (1 .- sum(s .* λvV[jprime] .* (Sxy[:, jprime] .> Sxy) for jprime in 1:length(v_y))) +
        sum(s .* (h_plus[:, jprime] * λvV') .* (Sxy .> Sxy[:, jprime]) for jprime in 1:length(v_y)) + 
       ( u_plus * λvV') .* (Sxy .≥ 0)
    )

    u_next = u_plus .* (1 .- sum(λvV[j] .* (Sxy[:, j] .≥ 0) for j in 1:length(v_y)))

    next::NamedTuple = (u = u_next, h = h_next)
    
    return next
    
end

export uh_next


## Need to get in initial conditions.
function uh_dynamics_sim(prod::Vector, u0::Vector, h0::Matrix, S::Array, l::Vector; Z = Z, X = X, Y = Y, δ = δ, s = s, α = α, c0 = c0, c1 = c1, ω = ω, p = p)
    T = length(prod)
    
    prod_zxy = valadd(Z, X, Y; p = p)
    GDP_zz = zeros(length(Z))
    next_zz = Vector{NamedTuple{}}(undef, length(Z))
    statet = Vector{Integer}(undef, T)
    GDP_zt = zeros(T)

    uxt = zeros(T, length(X))
    hxyt = zeros(T, length(X), length(Y))

    uxt[1, :] = u0
    hxyt[1, :, :] = h0
    statet[1] = Integer((length(Z) + 1)/2)

    for t in 2:T
        
        for zz in 1:length(Z)
            Sxy = S[zz, :, :]

            uplus = u_plus(uxt[t-1, :], hxyt[t-1, :, :], Sxy; δ = δ)
            hplus = h_plus(hxyt[t-1, :, :], Sxy; δ = δ)
            L = LT(uplus, hplus; s = s)
            Jy = J(uplus, hplus, L, Sxy; s = s)
            θ = tightness(Jy, L; α = α, c0 = c0, c1 = c1, ω = ω)
            v_y = vacy(Jy, θ; α = α, c0 = c0, c1 = c1, ω = ω)
    
            next_zz[zz] = uh_next(uplus, hplus, v_y, L, Sxy, l; α = α, ω = ω, s = s) 

            GDP_zz[zz] = sum(prod_zxy[zz, :, :] .* next_zz[zz].h)/(1 - sum(next_zz[zz].u))

        end

        statet[t] = argmin(abs.(13 .* GDP_zz .- prod[t]))
        next = next_zz[statet[t]]
        GDP_zt[t] = GDP_zz[statet[t]]

        uxt[t, :] = next.u
        hxyt[t, :, :] = next.h

    end

    return Dict(:uxt => uxt, :hxyt => hxyt, :statet => statet, :GDP_zt => GDP_zt)

end

export uh_dynamics_sim