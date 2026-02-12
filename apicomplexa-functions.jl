include("clustering.jl")
# %% Apicomplexa
leaq(a,b; kwargs...) = (a <= b) || isapprox(a, b; kwargs...)

function max_attained_at_least_twice(a,b,c)
    ma = maximum((a,b,c))
    mi = minimum((a,b,c))
    mid = a + b + c - ma - mi
    isapprox(ma, mid, atol=1e-12, rtol=sqrt(eps(mid)))
end

max_attained_precisely_twice(a,b,c) = max_attained_at_least_twice(a,b,c) && !(isapprox(a,b) && isapprox(a,c) && isapprox(b,c))

function is_ultrametric(m)
    n = size(m, 1)
    if !is_symmetric(m)
        return false
    end 
    for i in 1:n, j in i+1:n, k in j+1:n
        if !max_attained_at_least_twice(m[i,j], m[i,k], m[j,k])
            # println("Not ultrametric: d($i, $j) = $(m[i,j]), d($i, $k) = $(m[i,k]), d($j, $k) = $(m[j,k])")
            return false
        end
    end
    return true
end

function is_binary_tree(m)
    n = size(m, 1)
    if !is_symmetric(m)
        return false
    end 
    for i in 1:n, j in i+1:n, k in j+1:n
        if !max_attained_precisely_twice(m[i,j], m[i,k], m[j,k])
            # println("Not ultrametric: d($i, $j) = $(m[i,j]), d($i, $k) = $(m[i,k]), d($j, $k) = $(m[j,k])")
            return false
        end
    end
    return true
end

is_binary_tree(t::PhylogeneticTree) = is_binary_tree(cophenetic_matrix(t))

""" Make a phylogenetic tree equidistant by adding sufficient lengths to the edges
    adjacent to the leaves. All lengths of interiour edges remain the same. """
function make_equidistant(tree::PhylogeneticTree{T}) where T
    graph = adjacency_tree(tree)
    edge_lengths = tree.pm_ptree.EDGE_LENGTHS;
    height = tree.pm_ptree.NODE_HEIGHTS[1]
    leaf_indices = Dict(tree.pm_ptree.LEAVES[t]+1 => i for (i, t) in enumerate(taxa(tree)))
    m = cophenetic_matrix(tree)
    # The following function traverses the tree in DFS order, starting at the root (node 1).
    # It keeps track of the height `h` of the current node, and whenever it reaches a leaf,
    # it adds `height - h` to the corresponding row and column of `m`, extending the edge
    # to this leave by this difference.
    function f(v, h)
        for w in outneighbors(graph, v)
            h′ = h + edge_lengths[Edge(v, w)]
            if Oscar.outdegree(graph, w) == 0
                j = leaf_indices[w]
                l = convert(T, height - h′)
                for i in 1:size(m, 1)
                    m[i,j] += l
                    m[j,i] += l
                end
                m[j,j]  -= 2l
            else
                f(w, h′)
            end
        end
    end
    f(1, zero(height))
    @assert is_ultrametric(m) "The new distance matrix is not ultrametric."
    new_tree = phylogenetic_tree(m, taxa(tree))
    @assert is_equidistant(new_tree) "The new tree is not equidistant."
    return new_tree
end

""" Shift a cophenetic matrix onto the hyperplane H = {x | ∑ᵢxᵢ=0}. """
normalize_cophenetic_matrix(t::PhylogeneticTree{T}) where T = phylogenetic_tree(vech_to_matrix(vech(t) .- mean(vech(t))), taxa(t))
normalize_cophenetic_matrix(m::Union{QQMatrix, Matrix{Float64}}) = vech_to_matrix(vech(m) .- mean(vech(m)))
denormalize(t) = phylogenetic_tree(vech_to_matrix(1 .+ vech(t) .- minimum(vech(t))), taxa(t))
is_normalized(t) = isapproxzero(mean(vech(t)))
isapproxzero(t) = isapprox(t, zero(t); atol=1e-10)

# Alternative implementations of tropical median consensus tree, using `troipical_median` directly
include("make_tree_like_again.jl")
function tropical_median_consensus2(trees::AbstractVector{PhylogeneticTree{T}}) where T
    @assert all(is_equidistant.(trees)) "All input trees must be equidistant."
    t = only(unique(taxa.(trees)))
    mat = collect(transpose(stack(vech.(trees))))
    # mat .-= mean(mat, dims=2)
    sol = collect(convert(T, c) for c in Polymake.tropical.tropical_median(mat))::Vector{T}
    # sol .-= minimum(sol) # should probably do this (moving solution away from H), because Polymake might assume nonnegative entries for cophenetric matrix
    mat_sol = vech_to_matrix(sol)
    # @assert isapprox(sum(vech(mat))) "The resulting cophenetic matrix $mat does not lie on the hyperplane"
    if !is_ultrametric(mat_sol)
        println("ℹ️ The resulting cophenetic matrix $(vech(mat_sol)) is not ultrametric; try to make it tree-like again.")
        mat_sol = make_tree_like_again(mat_sol)
    end
    if !is_ultrametric(mat_sol)
        throw(mat_sol)
    end
    phylogenetic_tree(mat_sol, t)
end

# Alternative implementation of tropical median consensus tree, via linear programming
function tropical_median_consensus3(trees::AbstractVector{PhylogeneticTree{T}}) where T
    t = only(unique(taxa.(trees)))
    # Build and solve the corresponding LP directly, without using `tropical_median_consensus`.
    V = transpose(stack(vech.(trees)))
    m, n = size(V)
    V .-= mean(V, dims=2)
    M = T == Float64 ? zeros(m*n+2, m+n) : zero_matrix(QQ, m*n+2, m+n)
    for i in 1:m, j in 1:n
        M[(i-1) * n + j, i] = 1
        M[(i-1) * n + j, m + j] = 1
    end
    M[end-1, m+1:end] .=  1
    M[end  , m+1:end] .= -1
    v = [reshape(transpose(V), m*n); 0; 0]
    E = polyhedron(-M, -v)
    l = vcat(fill(n*one(T), m), fill(zero(T), n))
    @assert is_feasible(E)
    LP = linear_program(E, l; convention=:min)
    _, tx = solve_lp(LP)
    _ = tx[1:m]
    x = tx[m+1:end]
    mat = vech_to_matrix(x)
    # @assert iszero(sum(vech(mat))) "The resulting cophenetic matrix $mat does not lie on the hyperplane"
    # @assert is_ultrametric(mat) "The resulting cophenetic matrix $mat is not ultrametric."
    tree = phylogenetic_tree(mat, t)
    return tree
end

function height(t)
    @assert is_equidistant(t)
    t.pm_ptree.NODE_HEIGHTS[1]
end

function max_inner_height(t)
    @assert is_equidistant(t)
    maximum(t.pm_ptree.NODE_HEIGHTS[i] for i in 2:t.pm_ptree.N_NODES)
end

min_inner_depth(t) = height(t) - max_inner_height(t)

function min_inner_height(t)
    @assert is_equidistant(t)
    minimum(t.pm_ptree.NODE_HEIGHTS[i] for i in 2:t.pm_ptree.N_NODES if t.pm_ptree.NODE_DEGREES[i] != 1)
end

max_inner_depth(t) = height(t) - min_inner_height(t)

import Base: -, +
(-)(t::PhylogeneticTree, v) = phylogenetic_tree(vech_to_matrix(vech(t) .- 2*v), taxa(t))
(+)(t::PhylogeneticTree, v) = phylogenetic_tree(vech_to_matrix(vech(t) .+ 2*v), taxa(t))

function coarse_type(δ, t)
    n = length(taxa(t))
    m = cophenetic_matrix(t)

    # all (upper triangle) indices of m, ordered by entry
    u = [(i,j) for i in 1:n for j in i+1:n if m[i,j]/2 <= height(t) - δ]
    sort!(u, by=x->m[x[1], x[2]])

    # union-find data structure.
    repr = collect(1:n)
    represents = map(i->[i], 1:n)

    # invariant: all nodes with distance < m[i,j] belong to the same cluster.
    for (i,j) in u
        # merge the clusters of i and j
        if repr[i] != repr[j]
            append!(represents[repr[i]], represents[repr[j]])
            empty!(represents[repr[j]])
            repr[j] = repr[i]
        end
    end
    return Set(map(Set, filter(s->!isempty(s), represents)))
end

coarse_type(t) = coarse_type(min_inner_depth(t), t)

using Phylo
plot_phylo(t) = Plots.plot(parsenewick("($(newick(denormalize(t))[1:end-1]));")) 