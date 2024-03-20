function SurplusVFI(Z, X, Y, Π; MaxIter = 2000, tol = 10e-8, δ = 0.013, r = 0.05, p = [0.003, 2.053, -0.140, 8.035, -1.907, 6.596])
    
    s = flow_surplus(Z, X, Y; p)
    S = ((1 + r)/(r + δ)) .* s

    for iter in 1:max_iter
        S_next = [s[w, i, j] + ((1 - δ)/(1 + r)) + sum(Π[w, w_prime] * max(S[w_prime, i, j], 0) for w_prime in 1:length(Z)) for w in 1:length(Z), i in 1:length(X), j in 1:length(Y)]
        e = maximum(abs(S_next - S))
        if maximum(abs, S_next - S) < tol
            @info "The surplus function has succesfully converged after $iter iterations."
            return S_next
        elseif iter == max_iter
            @warn "After $iter iterations, the surplus function error remains larger than $tol."
            return S_next
        end
        S = S_next
    end



end
