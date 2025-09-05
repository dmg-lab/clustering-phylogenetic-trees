include("clustering.jl")
# %% Apicomplexa
leaq(a,b; kwargs...) = (a <= b) || isapprox(a, b; kwargs...)

is_ultrametric(m::Matrix{Float64}) = all(leaq(m[i,j], max(m[i,k], m[j,k]); atol=1e-10) for i in 1:size(m,1), j in 1:size(m,1), k in 1:size(m,1))

is_ultrametric(m::QQMatrix) = is_symmetric(m) && all(<=(m[i,j], max(m[i,k], m[j,k])) for i in 1:size(m,1), j in 1:size(m,1), k in 1:size(m,1))

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
            if outdegree(graph, w) == 0
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
normalize_cophenetic_matrix(m) = vech_to_matrix(vech(m) .- mean(vech(m)))

# Apicomplexa data set
samples = (open("R-Data/apicomplexa.txt")
     |> readlines
     |> Base.Fix2(getindex, load("R-Data/apicomplexa_path_subset.txt")) # Subset that causes the problem to appear
    .|> (s -> phylogenetic_tree(Float64, s))
    .|> make_equidistant
)

# reduce sample size for testing
samples = samples[1:25]

# Tropical median consensus tree via Andrei's implementation
@time mt1 = tropical_median_consensus(samples)
@assert is_equidistant(mt1)
nm1 = normalize_cophenetic_matrix(cophenetic_matrix(mt1))
d1 = Float64(sum(d.(samples, Ref(mt1))))

# Alternative implementations of tropical median consensus tree, using `troipical_median` directly
function tropical_median_consensus2(trees::AbstractVector{PhylogeneticTree{T}}) where T
    t = only(unique(taxa.(trees)))
    mat = collect(transpose(stack(vech.(trees))))
    mat .-= mean(mat, dims=2)
    sol = collect(convert(T, c) for c in Polymake.tropical.tropical_median(mat))::Vector{T}
    sol .-= minimum(sol)
    mat = vech_to_matrix(sol)
    @assert is_ultrametric(mat) "The resulting cophenetic matrix is not ultrametric: $mat"
    phylogenetic_tree(mat, t)
end

@time mt2 = tropical_median_consensus2(samples)
@assert is_equidistant(mt2)
nm2 = normalize_cophenetic_matrix(cophenetic_matrix(mt2))
d2 = Float64(sum(d.(samples, Ref(mt2))))

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
    tree = phylogenetic_tree(vech_to_matrix(x), t)
    return tree
end

@time mt3 = tropical_median_consensus3(samples)
nm3 = normalize_cophenetic_matrix(cophenetic_matrix(mt3))
d3 = sum(d.(samples, Ref(mt3)))

isapprox.(nm1, nm2)
isapprox.(nm1, nm3)
isapprox.(nm2, nm3)
isapprox(d1, d2)
isapprox(d1, d3)
isapprox(d2, d3)

# More debugging stuff (Lena)
# ===========================
# reduce sample size
for i in 1:84-27
    a, b = i, 27+i
    tmc = tropical_median_consensus(s[a:b])

    V = transpose(stack(vech.(s[a:b])))
    m, n = size(V)
    M = zeros(m*n+2, m+n)
    for i in 1:m, j in 1:n
        M[(i-1) * n + j, i] = 1.0
        M[(i-1) * n + j, m + j] = 1.0
    end
    M[end-1, m+1:end] .=  1
    M[end  , m+1:end] .= -1
    v = [reshape(transpose(V), m*n); 0; 0]
    E = polyhedron(-M, -v)
    l = vcat(n*ones(m), zeros(n))
    LP = linear_program(E, l; convention=:min)
    r, tx = solve_lp(LP)
    t = tx[1:m]
    x = tx[m+1:end]
    tree = phylogenetic_tree(vech_to_matrix(x), taxa(s[1]))


    if sum(d.(s[a:b], Ref(tmc))) - sum(d.(s[a:b], Ref(tree))) > 1
        println(i)
    end
end

# reconstruct polymake side for a,b = 1,28
# Note that if we compare the resulting cophenetic matrices
cophenetic_matrix(tree)
cophenetic_matrix(tmc)
# particularly the entries [6,8] and [7,8] in the c. matrix of tmc are significantly bigger than those of tree (even after subtracting the constant factor that Andrei adds)
S = zeros(28,28) # essentially V
for k in 1:8
    for (i,j) in combinations(1:8,2)
            S[k,8*(i-1) - Int(i*(i+1)//2)+ j] = cophenetic_matrix(s[k])[i,j]
    end
end
# As expected, when comparing the result of the tropical_median function with x above the last two entries diverge significantly after making up for the
tm = Vector(Polymake.call_function(:tropical,:tropical_median, S))
filter(i -> abs(tm[i] - x[i]) > 1, 1:28)

# debugging tropical_median
# 1. computing tropical vertices
m, n = 28, 28
supply = n*ones(m); demand = m*ones(n)
flowMatrix = Polymake.call_function(:graph,:optimal_transport_plan, -S, supply, demand)
trop_vert = Polymake.call_function(:tropical, :facets_matrix, S, flowMatrix)
# 2. computing the average of the tropical vertices
tm = zeros(n)
r, c = size(trop_vert)
for i in 1:c
    tm[i] = 1/r*(sum(trop_vert[:,i]))
end
# The result is similar to the result of the tropical_median function in polymake (up to adding some multiple of the all ones vector)
# The problem must be in the computation of the tropical vertices then