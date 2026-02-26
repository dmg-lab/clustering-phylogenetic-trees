"""
    repeat_until_success(f, args...; kwargs...)

Runs `f(args...; kwargs...)` until no exception is raised.
Of course, this only makes sense for non-deterministic functions.
"""
function repeat_until_success(f, args...; kwargs...)
    while true
        try
            return f(args...; kwargs...)
        catch e
            println("⚠️ Failed; try again")
            continue
        end
    end
end

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

"""
    tropical_median_consensus_div(samples)

If `length(samples)` is divisible by the dimension,
then (for numeric instability reasons), one sample is excluded,
before `invoking tropical_median_consensus2`.
"""
function tropical_median_consensus_div(samples)
    n = length(taxa(samples[1]))
    if length(samples) % (n*(n-1)//2) == 0
        println("⚠️ Excluding one sample for divisibility reasons.")
        r = Int(rand(UInt) % length(samples)) + 1
        tropical_median_consensus2(samples[[1:r; r+2:end]])
    else
        tropical_median_consensus2(samples[1:end])
    end
end
