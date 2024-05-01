using thesis_code, Plots

grid = grids()
X = grid[:X]
Z = grid[:Z]

function heat_sub(sub::Vector)
    sub_array = [max((sub[1] + sub[2] * Z[w] + sub[3] * X[i] + sub[4] * Z[w] * X[i]), 0) for i in 1:lastindex(X), w in 1:lastindex(Z)]
    p = heatmap(X, Z, sub_array)
    
    return p
end

sub = parse.(Float64, readlines("output/optsub.txt"))

head_sub(sub)

sub_ui = parse.(Float64, readlines("output/optsub_ui.txt"))

head_sub(sub_ui)