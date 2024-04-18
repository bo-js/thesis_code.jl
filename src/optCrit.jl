function optCrit(g, τ, grid::Dict, params::Dict; T = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52)
    
    # Get Params
    p = params[:p]

    # Define Grids 
    X = grid[:X]
    Y = grid[:Y]
    Z = grid[:Z]
    

    va = valadd(Z, X, Y; p = p)
    b = homeprod(X, Y; p = p)

    sub = [((Z[w] < 1)*(g[1] * Y[j] + g[2] * Z[w] + g[3] * Y[j] * Z[w] - τ)) * va[w, i, j] for w in 1:lastindex(Z), i in 1:lastindex(X), j in 1:lastindex(Y)]    

    fs = flow_surplus(Z, X, Y, sub; p = p)

    sim = sim_draws(grid, params, fs; T = T, burn = burn, draw = draw, r = r )
    
    crit = sum((sum(sim.hxyt[t, :, :] .* va[sim.statet[t], :, :]) + dot(sim.uxt[t, :], b)) for t in burn+1:T+burn)
    
    deficit = sum(sum(sim.hxyt[t, :, :] * sub[sim.statet[t], :, :]) for t in burn+1:T+burn)

    return (c = crit, d = deficit)
end

export optCrit