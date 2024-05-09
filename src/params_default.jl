function params_default(; errata = false)
    
    if errata == false
        d = (
            ## Matching
        α = 0.497,
        ω = 1/2,
        ## Search Intensity
        s = 0.027,
        ## Vacancy Costs
        c0 = 0.028,
        c1 = 0.084,
        ## Exogeneous Separation
        δ = 0.013,
        ## Productivity Shocks
        σ = 0.071,
        ρ = 0.999,
        ## Worker Heterogeneity
        β1 = 2.148,
        β2 = 12.001,
        ## Value Added
        p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596]
        )
    else
        d = (
            ## Matching
        α = 0.645,
        ω = 1/2,
        ## Search Intensity
        s = 0.023,
        ## Vacancy Costs
        c0 = 0.049,
        c1 = 0.007,
        ## Exogeneous Separation
        δ = 0.013,
        ## Productivity Shocks
        σ = 0.072,
        ρ = 0.999,
        ## Worker Heterogeneity
        β1 = 2.329,
        β2 = 17.417,
        ## Value Added
        p = [0.000, 0.991, -0.126, 9.042, -0.425, 3.391]
        )
    end
    return d
end

export params_default