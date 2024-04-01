function u_plus(ux::Vector, hxy::Matrix, Sxy::Matrix; δ::Number = δ)
   return ux .+ sum(((Sxy[:, j] .< 0) + δ.*(Sxy[:, j] .≥ 0)) .* hxy[:, j] for j in 1:length(Sxy[1, :])) 
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

function uh_next(u_plus::Vector, h_plus::Matrix, v_y::Vector, L::Number, Sxy::Matrix; α::Number = α, ω::Number = ω, s::Number = s, calc_u = true::Bool)

    V = sum(v_y)
    λ = min(α * L^ω * V^(1 - ω), L, V)/L
    λvV = λ .* v_y ./V

    if calc_u == true
        u_next = max.(u_plus .* (1 .- sum(λvV[j] .* (Sxy[:, j].≥0) for j in 1:length(v_y))), 0)
    end

    h_next = max.(
        h_plus .* (1 .- sum((s * λvV[jprime]) .* (Sxy[:, jprime] .> Sxy) for jprime in 1:length(v_y))) +
        sum(h_plus[:, jprime] .* s .* λvV .* (Sxy .> Sxy[:, jprime]) for jprime in 1:length(v_y)) +
        (u_plus * λvV') .* (Sxy .≥ 0), 0
    )

    if calc_u == true
        
        next::NamedTuple = (u = u_next, h = h_next)
    
        return next
    else
        return h_next
    end
end

export uh_next


## Need to get in initial conditions.
function uh_dynamics_sim(prod::Vector, u0::Vector, h0::Matrix, S::Array; Z = Z, X = X, Y = Y)
    T = length(prod)
    
    prod_zxy = valadd(Z, X, Y; p = p)
    GDP_z = zeros(Nz)
    uh_zz = zeros(Nz)
    statet = zeros(T)
    GDP_zt = zeros(T)

    uxt = zeros(T, Nx)
    hxyt = zeros(T, Nx, Ny)

    uxt[1, :] = u0
    hxyt[1, :, :] = h0

    for t in 2:T
        
        for zz in 1:Nz
            uh_zz[zz] = uh_next(uxt, hxyt, S, zz, t - 1)
            
            GDP_z[zz] = sum(prod_zxy[zz, :, :] .* uh_zz[zz].h_next)/(1 - sum(uh_zz[zz].u_next))

        end

        statet[t] = argmin(abs.(13 .* GDP_z .- prod[t]))
        next = uh_zz[statet[t]]
        GDP_zt[t] = GDP_z[state[t]]

        uxt[1, :] = next.u_next
        hxyt[1, :, :] = next.h_next

    end

    return Dict(uxt => :uxt, hxyt => :hxyt, statet => :statet, GDP_zt => :GDP_zt)
end

export uh_dynamics_sim