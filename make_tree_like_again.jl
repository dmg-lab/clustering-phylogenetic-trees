function make_tree_like_again(mat)
    m = copy(mat)
    # @assert size(m, 1) == size(m, 2) && is_symmetric(m)
    n = size(m, 1)
    u = [(i,j) for i in 1:n for j in i+1:n]
    sort!(u, by=x->m[x[1], x[2]])

    repr = collect(1:n)
    represents = map(i->[i], 1:n)

    for (i,j) in u
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
    if max_delta > 1e-10
        println("⚠️ The maximum change to make the matrix tree-like is $max_delta, which is larger than the tolerance of 1e-10.")
    end
    @assert is_ultrametric(m) "The resulting matrix is not ultrametric: $m"
    return m
end
