function params_default()
    d = Dict(
        ## Matching
    :α => 0.497,
    :ω => 1/2,
    ## Search Intensity
    :s => 0.027,
    ## Vacancy Costs
    :c0 => 0.028,
    :c1 => 0.084,
    ## Exogeneous Separation
    :δ => 0.013,
    ## Productivity Shocks
    :σ => 0.071,
    :ρ => 0.999,
    ## Worker Heterogeneity
    :β1 => 2.148,
    :β2 => 12.001,
    ## Value Added
    :p => [0.003, 2.053, -0.140, 8.035, -1.907, 6.596]
    )
    return d
end

export params_default