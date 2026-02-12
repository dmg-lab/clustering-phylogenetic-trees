"""
    make_tree_like_again(m)

Given a symmetric matrix representing a finite matrix m,
return a matrix m' such that
a) m' .<= m
b) maximum(m'[i,j], m'[i,k], m'[j,k]) is attained twice.
This is done by setting the maximum to the middle term.
"""
function make_tree_like_again(mat)
    @assert size(mat, 1) == size(mat, 2) && is_symmetric(mat)
    m = deepcopy(mat)
    n = size(m, 1)

    # all (upper triangle) indices of m, ordered by entry
    u = [(i,j) for i in 1:n for j in i+1:n]
    sort!(u, by=x->m[x[1], x[2]])

    # union-find data structure.
    repr = collect(1:n)
    represents = map(i->[i], 1:n)

    # invariant: all nodes with distance < m[i,j] belong to the same cluster.
    for (i,j) in u
        # merge the clusters of i and j
        if repr[i] != repr[j]
            for k in represents[repr[i]], l in represents[repr[j]]
                m[k,l] = m[l,k] = m[i,j]
            end
            append!(represents[repr[i]], represents[repr[j]])
            empty!(represents[repr[j]])
            repr[j] = repr[i]
        end
    end
    max_delta = maximum(abs.(m - mat))
    if max_delta > 1e-7
        println("⚠️ The maximum change to make the matrix tree-like is $max_delta, which is larger than the tolerance of 1e-10.")
    end
    @assert is_ultrametric(m) "The resulting matrix is not ultrametric: $m"
    return m
end

make_tree_like_again(t::PhylogeneticTree) = phylogenetic_tree(make_tree_like_again(cophenetic_matrix(t)), taxa(t))
