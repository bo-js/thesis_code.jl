function ES(S, Π)
    ExpS = zeros(size(S))
    for j in 1:length(S[1, 1, :])
        ExpS[:, :, j] = Π * max.(S[:, :, j], 0)
    end
    return ExpS
end

function SurplusVFI(Z, X, Y, Π; MaxIter = 2000, tol = 10e-8, δ = 0.013, r = 0.05, p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    
    s = flow_surplus(Z, X, Y; p)
    S = ((1 + r)/(r + δ)) .* s

    for iter in 1:MaxIter
        S_next = s + ((1 - δ)/(1 + r)) * ES(S, Π)
        
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
