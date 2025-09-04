using Oscar, Combinatorics

include("clustering.jl")

# %% Consensus tree bug
f(m, t) = phylogenetic_tree(m, t)
f(m, t) = phylogenetic_tree(matrix(QQ, convert(Matrix{Rational},m)), t)
f(m, t) = make_equidistant(phylogenetic_tree(matrix(QQ, convert(Matrix{Rational},m)), t))
t, s = load("consensus_tree_bug.json")
t, s = f(t...), [f(a, b) for (a, b) in s]
tmc = tropical_median_consensus(s)
sum(d.(s, Ref(tmc)))

# %% Reduce sample size & map to H
V = transpose(stack(vech.(s[1:28])))
m, n = size(V)
hDiff = sum.([V[i,:]/n for i in 1:m])
for i in 1:m
    for j in 1:n
        V[i,j] -= hDiff[i]
    end
end
V = Matrix(V)

# Find tropical median consensus of samples after there where mapped to H
sample_on_H = [phylogenetic_tree(vech_to_matrix(V[i,:]),taxa(s[i])) for i in 1:nrows(V)]
make_equidistant.(sample_on_H)

flowMatrix = Polymake.call_function(:graph,:optimal_transport_plan, -V, supply, demand)
V_facets = facets_matrix(V, Matrix(flowMatrix))
sum_avg = 1/ncols(facets)*sum.([facets[i,:] for i in 1:nrows(facets)])
for i in 1:nrows(V)
    V[i,:] -= sum_avg[i]*ones(ncols(V))
end
1/28*sum([V[i,:] for i in 1:nrows(V)])