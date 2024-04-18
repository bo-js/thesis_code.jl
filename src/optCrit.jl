function optCrit(g, τ; grid::Dict = grids(), params::Dict = params_default(), T = 5000, burn = 1000, draw = rand(burn+T, 1),  r = 0.05/52, unconstr = false)
    
    if unconstr == true
        if τ > 1
            return 0
        elseif τ < 0
            return 0
        end
    end
    
    # Get Params
    p = params[:p]

    # Define Grids 
    X = grid[:X]
    Y = grid[:Y]
    Z = grid[:Z]
    

    va = valadd(Z, X, Y; p = p)
    b = homeprod(X, Y; p = p)

    sub = repeat((Z .< 1), 1, length(X), length(Y)) .* (g .* va.^(-1)) .- τ .* va

    fs = flow_surplus(va+sub, b)

    sim = sim_draws(grid, params, fs; T = T, burn = burn, draw = draw, r = r )
    
    crit = sum((sum(sim.hxyt[t, :, :] .* (va[sim.statet[t], :, :] + sub[sim.statet[t], :, :])) + dot(sim.uxt[t, :], b)) for t in burn+1:T+burn)
    
    deficit = sum(sum(sim.hxyt[t, :, :] * sub[sim.statet[t], :, :]) for t in burn+1:T+burn)

    if unconstr == false
        return (c = crit, d = deficit)
    elseif deficit ≤ 0
        return crit
    else
        return 0
    end
end

export optCrit