function ES(S::Array, Π::Matrix; threshold = 0)
    ExpS = zeros(size(S))
    for j in 1:length(S[1, 1, :])
        ExpS[:, :, j] = Π * max.(S[:, :, j], threshold)
    end
    return ExpS
end

function ES(S::Array, Π::Matrix, thresholds::Array)
    ExpS = zeros(size(S))
    for j in 1:length(S[1, 1, :])
        ExpS[:, :, j] = Π * max.(S[:, :, j], thresholds[:, :, j])
    end
    return ExpS
end

function SurplusVFI(Z::Vector, X::Vector, Y::Vector, Π::Matrix; threshold = 0, MaxIter = 2000, tol = 10e-8, δ = δ, r = r, p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    
    s = flow_surplus(Z, X, Y, p = p)
    S = ((1 + r)/(r + δ)) .* s

    for iter in 1:MaxIter
        S_next = s + ((1 - δ)/(1 + r)) * ES(S, Π; threshold = threshold)
        
        if maximum(abs, S_next - S) < tol
            @info "The surplus function has succesfully converged after $iter iterations."
            return S_next
        elseif iter == MaxIter
            @warn "After $iter iterations, the surplus function error remains larger than $tol."
            return S_next
        end
        S = S_next
    end

end

function SurplusVFI(flowsurp::Array, Π;threshold = 0, MaxIter = 5000, tol = 10e-8, δ = δ, r = r, p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    
    S = ((1 + r)/(r + δ)) .* flowsurp

    for iter in 1:MaxIter
        S_next = flowsurp + ((1 - δ)/(1 + r)) * ES(S, Π; threshold = threshold)
        
        if maximum(abs, S_next - S) < tol
            @info "The surplus function has succesfully converged after $iter iterations."
            return S_next
        elseif iter == MaxIter
            @warn "After $iter iterations, the surplus function error remains larger than $tol."
            return S_next
        end
        S = S_next
    end

end

function SurplusVFI(flowsurp::Array, Π, thresholds; MaxIter = 5000, tol = 10e-8, δ = δ, r = r, p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    
    S = ((1 + r)/(r + δ)) .* flowsurp

    for iter in 1:MaxIter
        S_next = flowsurp + ((1 - δ)/(1 + r)) * ES(S, Π, thresholds)
        
        if maximum(abs, S_next - S) < tol
            @info "The surplus function has succesfully converged after $iter iterations."
            return S_next
        elseif iter == MaxIter
            @warn "After $iter iterations, the surplus function error remains larger than $tol."
            return S_next
        end
        S = S_next
    end

end

function H0iter(Z, X, Y, S, l; MaxIter = 5000, tol = 10e-5, δ = δ, s = s, α = α, c0 = c0, c1 = c1, ω = ω)
    z = Integer((length(Z)+1)/2)
    Sxy = S[z, :, :]

    ## Initialize workers as uniformly distributed over job types
    h0_xy =  [l[i]/length(Y) for i in 1:length(X), j in 1:length(Y)]
    u0_x = zeros(length(X))
    for iter in 1:MaxIter

        uplus = u_plus(u0_x, h0_xy, Sxy; δ = δ)
        hplus = h_plus(h0_xy, Sxy; δ = δ)
        L = LT(uplus, hplus; s = s)
        Jy = J(uplus, hplus, L, Sxy; s = s)
        θ = tightness(Jy, L; α = α, c0 = c0, c1 = c1, ω = ω)
        v_y = vacy(Jy, θ; α = α, c0 = c0, c1 = c1, ω = ω)

        next = uh_next(uplus, hplus, v_y, L, Sxy, l; α = α, ω = ω, s = s)

        if maximum(abs, next.h - h0_xy) < tol
            @info "H0 has succesfully converged after $iter iterations."
            return next
        elseif iter == MaxIter
            @warn "After $iter iterations, the error remains larger than $tol."
            return next
        end
        h0_xy = next.h
        u0_x = next.u
    end

end
