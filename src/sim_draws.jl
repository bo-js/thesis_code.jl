
function sim_draws(grd::Dict, params::NamedTuple; threshold = 0, T::Integer = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52)

    α = params.α
    ω = params.ω
    s = params.s
    c0 = params.c0
    c1 = params.c1
    δ = params.δ

   # Define Grids 
    X = grd[:X]
    Y = grd[:Y]
    Z = grd[:Z]
    Π = grd[:Π]
    l = grd[:l]

    # Find Surplus Function and Initial Conditions
    S = SurplusVFI(Z, X, Y, Π;threshold = threshold, δ = δ, r = r)
    init = H0iter(Z, X, Y, S, l; δ = δ, s=s, α = α, c0 = c0, c1 = c1, ω = ω)

    statet = Vector{Integer}(undef, T+burn)
    uxt = zeros(T+burn, length(X))
    hxyt = zeros(T+burn, length(X), length(Y))
    vacyt = zeros(T+burn, length(Y))


    uxt[1, :] = init.u
    hxyt[1, :, :] = init.h
    i = Integer(round(length(Z)/2))
    statet[1] = i

    for t in 2:(T+burn)
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), length(Z))
        statet[t] = i

        Sxy = S[i, :, :]

        uplus = u_plus(uxt[t-1, :], hxyt[t-1, :, :], Sxy; δ = δ)
        hplus = h_plus(hxyt[t-1, :, :], Sxy; δ = δ)
        L = LT(uplus, hplus; s = s)
        Jy = J(uplus, hplus, L, Sxy; s = s)
        θ = tightness(Jy, L; α = α, c0 = c0, c1 = c1, ω = ω)
        v_y = vacy(Jy, θ; α = α, c0 = c0, c1 = c1, ω = ω)   
        vacyt[t, :] = v_y
 
        
        next = uh_next(uplus, hplus, v_y, L, Sxy, l; α = α, ω = ω, s = s)
        uxt[t, :] = next.u
        hxyt[t, :, :] = next.h
    end
    
    return (statet = statet, uxt = uxt, hxyt = hxyt, vacyt = vacyt)

end

function sim_draws(grd::Dict, params::NamedTuple, flowsurp::Array; threshold = 0, T::Integer = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52)

    α = params.α
    ω = params.ω
    s = params.s
    c0 = params.c0
    c1 = params.c1
    δ = params.δ

   # Define Grids 
    X = grd[:X]
    Y = grd[:Y]
    Z = grd[:Z]
    Π = grd[:Π]
    l = grd[:l]

    # Find Surplus Function and Initial Conditions
    S = SurplusVFI(flowsurp, Π; δ = δ, r = r, threshold = threshold)

    init = H0iter(Z, X, Y, S, l; δ = δ, s=s, α = α, c0 = c0, c1 = c1, ω = ω)

    statet = Vector{Integer}(undef, T+burn)
    uxt = zeros(T+burn, length(X))
    hxyt = zeros(T+burn, length(X), length(Y))
    vacyt = zeros(T+burn, length(Y))

    uxt[1, :] = init.u
    hxyt[1, :, :] = init.h
    i = Integer(round(length(Z)/2))
    statet[1] = i

    for t in 2:(T+burn)
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), length(Z))
        statet[t] = i

        Sxy = S[i, :, :]

        uplus = u_plus(uxt[t-1, :], hxyt[t-1, :, :], Sxy; δ = δ)
        hplus = h_plus(hxyt[t-1, :, :], Sxy; δ = δ)
        L = LT(uplus, hplus; s = s)
        Jy = J(uplus, hplus, L, Sxy; s = s)
        θ = tightness(Jy, L; α = α, c0 = c0, c1 = c1, ω = ω)
        v_y = vacy(Jy, θ; α = α, c0 = c0, c1 = c1, ω = ω)    
        vacyt[t, :] = v_y

        next = uh_next(uplus, hplus, v_y, L, Sxy, l; α = α, ω = ω, s = s)
        uxt[t, :] = next.u
        hxyt[t, :, :] = next.h
    end
    
    return (statet = statet, uxt = uxt, hxyt = hxyt, vacyt = vacyt)

end

function sim_draws(grd::Dict, params::NamedTuple, flowsurp::Array, init::NamedTuple; threshold = 0, T::Integer = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52)

    α = params.α
    ω = params.ω
    s = params.s
    c0 = params.c0
    c1 = params.c1
    δ = params.δ

   # Define Grids 
    X = grd[:X]
    Y = grd[:Y]
    Z = grd[:Z]
    Π = grd[:Π]
    l = grd[:l]

    # Find Surplus Function and Initial Conditions
    S = SurplusVFI(flowsurp, Π; threshold = threshold, δ = δ, r = r)

    statet = Vector{Integer}(undef, T+burn)
    uxt = zeros(T+burn, length(X))
    hxyt = zeros(T+burn, length(X), length(Y))
    vacyt = zeros(T+burn, length(Y))

    uxt[1, :] = init.u
    hxyt[1, :, :] = init.h
    i = Integer(round(length(Z)/2))
    statet[1] = i

    for t in 2:(T+burn)
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), length(Z))
        statet[t] = i

        Sxy = S[i, :, :]

        uplus = u_plus(uxt[t-1, :], hxyt[t-1, :, :], Sxy; δ = δ)
        hplus = h_plus(hxyt[t-1, :, :], Sxy; δ = δ)
        L = LT(uplus, hplus; s = s)
        Jy = J(uplus, hplus, L, Sxy; s = s)
        θ = tightness(Jy, L; α = α, c0 = c0, c1 = c1, ω = ω)
        v_y = vacy(Jy, θ; α = α, c0 = c0, c1 = c1, ω = ω)    
        vacyt[t, :] = v_y

        next = uh_next(uplus, hplus, v_y, L, Sxy, l; α = α, ω = ω, s = s)
        uxt[t, :] = next.u
        hxyt[t, :, :] = next.h
    end
    
    return (statet = statet, uxt = uxt, hxyt = hxyt, vacyt = vacyt)

end

function sim_draws(grd::Dict, params::NamedTuple, flowsurp::Array, init::NamedTuple, thresholds::Array; T::Integer = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52)

    α = params.α
    ω = params.ω
    s = params.s
    c0 = params.c0
    c1 = params.c1
    δ = params.δ

   # Define Grids 
    X = grd[:X]
    Y = grd[:Y]
    Z = grd[:Z]
    Π = grd[:Π]
    l = grd[:l]

    # Find Surplus Function and Initial Conditions
    S = SurplusVFI(flowsurp, Π, thresholds; δ = δ, r = r)

    statet = Vector{Integer}(undef, T+burn)
    uxt = zeros(T+burn, length(X))
    hxyt = zeros(T+burn, length(X), length(Y))
    vacyt = zeros(T+burn, length(Y))

    uxt[1, :] = init.u
    hxyt[1, :, :] = init.h
    i = Integer(round(length(Z)/2))
    statet[1] = i

    for t in 2:(T+burn)
        i = min(1 + sum(draw[t] .> cumsum(Π[i, :])), length(Z))
        statet[t] = i

        Sxy = S[i, :, :]

        uplus = u_plus(uxt[t-1, :], hxyt[t-1, :, :], Sxy; δ = δ)
        hplus = h_plus(hxyt[t-1, :, :], Sxy; δ = δ)
        L = LT(uplus, hplus; s = s)
        Jy = J(uplus, hplus, L, Sxy; s = s)
        θ = tightness(Jy, L; α = α, c0 = c0, c1 = c1, ω = ω)
        v_y = vacy(Jy, θ; α = α, c0 = c0, c1 = c1, ω = ω)    
        vacyt[t, :] = v_y

        next = uh_next(uplus, hplus, v_y, L, Sxy, l; α = α, ω = ω, s = s)
        uxt[t, :] = next.u
        hxyt[t, :, :] = next.h
    end
    
    return (statet = statet, uxt = uxt, hxyt = hxyt, vacyt = vacyt)

end


export sim_draws