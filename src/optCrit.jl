function optCrit(g1, g2, g3, g4, ui, τ; grid::Dict = grids(), params::NamedTuple = params_default(), T = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52, unconstr = false)
    
    if unconstr == true
        if τ > 1 
            return - Inf
        elseif τ < 0 || g2 < 0
            return - Inf
        end
    end
    
    # Get Params
    p = params.p
    c0 = params.c0
    c1 = params.c1

    # Define Grids 
    X = grid[:X]
    Y = grid[:Y]
    Z = grid[:Z]
    

    va = valadd(Z, X, Y; p = p)
    b = homeprod(X, Y; p = p)
    uix = uipayments(X, Y, ui; p = p)

    ##sub = repeat((Z .< 1), 1, length(X), length(Y)) .* ((g1 .- g2 .* va)) .- τ .* va
    sub = [max((g1 - g2 * Z[w] - g3 * X[i] - g4 * Z[w] * X[i]), 0) .- τ .* va[w, i, j] for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]

    fs = flow_surplus(va .+ sub, b .+ uix)

    sim = sim_draws(grid, params, fs; T = T, burn = burn, draw = draw, r = r )
    
    crit = (sum(sum(sim.hxyt[t, :, :] .* (va[sim.statet[t], :, :] + sub[sim.statet[t], :, :])) 
        + dot(sim.uxt[t, :], b + uix) - (c0/(1 + c1)) * sum(sim.vacyt[t, :] .^ (1 + c1)) for t in burn+1:T+burn))/T
    
    deficit = sum(sum(sim.hxyt[t, :, :] * sub[sim.statet[t], :, :]) + dot(sim.uxt[t, :], uix) for t in burn+1:T+burn)

    if unconstr == false
        return (c = crit, d = deficit)
    elseif deficit ≤ 0
        return crit
    else
        return - Inf
    end
end

function optCrit(g1, g2, g3, g4, ui, τ, init::NamedTuple; grid::Dict = grids(), params::NamedTuple = params_default(), T = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52, unconstr = false)
    
    if unconstr == true
        if τ > 1 
            return - Inf
        elseif τ < 0 || g2 < 0
            return - Inf
        end
    end
    
    # Get Params
    p = params.p
    c0 = params.c0
    c1 = params.c1

    # Define Grids 
    X = grid[:X]
    Y = grid[:Y]
    Z = grid[:Z]
    

    va = valadd(Z, X, Y; p = p)
    b = homeprod(X, Y; p = p)
    uix = uipayments(X, Y, ui; p = p)

    ##sub = repeat((Z .< 1), 1, length(X), length(Y)) .* ((g1 .- g2 .* va)) .- τ .* va
    sub = [max((g1 - g2 * Z[w] - g3 * X[i] - g4 * Z[w] * X[i]), 0) .- τ .* va[w, i, j] for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]

    fs = flow_surplus(va .+ sub, b .+ uix)

    sim = sim_draws(grid, params, fs, init; T = T, burn = burn, draw = draw, r = r )
    
    crit = (sum(sum(sim.hxyt[t, :, :] .* (va[sim.statet[t], :, :] + sub[sim.statet[t], :, :])) 
        + dot(sim.uxt[t, :], b + uix) - (c0/(1 + c1)) * sum(sim.vacyt[t, :] .^ (1 + c1)) for t in burn+1:T+burn))/T
    
    deficit = sum(sum(sim.hxyt[t, :, :] * sub[sim.statet[t], :, :]) + dot(sim.uxt[t, :], uix) for t in burn+1:T+burn)

    if unconstr == false
        return (c = crit, d = deficit)
    elseif deficit ≤ 0
        return crit
    else
        return - Inf
    end
end

export optCrit