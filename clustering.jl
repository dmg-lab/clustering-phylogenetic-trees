using Oscar, LinearAlgebra, Plots, JSON, OrderedCollections, Random, StatsBase, Printf
import Oscar: PhylogeneticTree
const col = palette(:tab10)
const col0 = palette([:black])

""" Like `first`, but returns `default` if the iterator is empty."""
function firstd(iter, default=nothing)
    for x in iter
        return x
    end
    return default
end

""" Convert a square matrix to the vector of its upper triangular part, in row-major order."""
function vech(A::AbstractMatrix{T}) where T
    m = LinearAlgebra.checksquare(A)
    v = Vector{T}(undef, m*(m - 1) ÷ 2)
    k = 0
    for i in 1:m, j in i+1:m
        @inbounds v[k += 1] = A[i,j]
    end
    return v
end

""" Convert a vector to a symmetric matrix, in row-major order, with diagonals set to zero."""
function vech_to_matrix(v::AbstractVector{T}) where T
    n = Int((sqrt(1 + 8*length(v)) + 1) ÷ 2)
    @assert n * (n - 1) ÷ 2 == length(v) "Vector length must be n(n-1)/2 for some n."
    A = zeros(T, (n, n))
    k = 0
    for i in 1:n, j in i+1:n
        @inbounds A[j, i] = A[i, j] = v[k += 1]
    end
    return A
end

""" Convert a PhylogeneticTree to a vector in ℝᵉ/𝟏ℝ, where e = choose(#leaves, 2)."""
vech(t::PhylogeneticTree) = vech(cophenetic_matrix(t))

""" Asymmetric tropical distance. """
function d(x::Vector, y::Vector)
    sum(y .- x) - length(x)*minimum(y .- x)
end

function d(x::PhylogeneticTree, y::PhylogeneticTree)
    nx = length(taxa(x))
    ny = length(taxa(y))
    @assert nx == ny "Cannot compare trees with different numbers of leaves."
    d(vech(x), vech(y))
end

function normalize(newick_string)
    i = 0
    r(_) = "t$(i += 1)"
    replace(newick_string, r"t[0-9]+" => r)
end

strip_parentheses(newick_string) = replace(newick_string, r"^\((.*)\);$" => s"\1;")

read_newick_json(filename, do_normalize=false) = (
        open(filename)
     |> readlines
     |> join
     |> JSON.parse
    .|> first
    .|> strip_parentheses
    .|> (do_normalize ? normalize : identity)
    .|> Base.Fix1(phylogenetic_tree, Float64)
)

samples = read_newick_json("R-Data/newick-100-4.json", true)

Random.seed!(1)

function farthest_point_sampling_rand(samples, k; d=d)
    cs = Vector{PhylogeneticTree}(undef, k)
    cs[1] = rand(samples)
    ds = d.(Ref(cs[1]), samples)
    for i in 2:k
        cs[i] = sample(samples, pweights(ds))
        ds = min.(ds, d.(Ref(cs[i]), samples))
    end
    return cs
end

function farthest_point_sampling_strict(samples, k; d=d)
    cs = Vector{PhylogeneticTree}(undef, k)
    cs[1] = rand(samples)
    ds = d.(Ref(cs[1]), samples)
    for i in 2:k
        cs[i] = samples[argmax(ds)]
        ds = min.(ds, d.(Ref(cs[i]), samples))
    end
    return cs
end

""" k-means-clustering w.r.t. the distance function `d`.
    Either provide `centroids` as a vector of trees, a vector of indices into `samples`, or an integer.
    In the latter case, the function samples `centrs` many centroids at random."""
function cluster(samples, centrs::Union{Int, Vector{Int}, Vector{PhylogeneticTree}}; d=d)
    centroids = centrs isa Int ? farthest_point_sampling_rand(samples, centrs; d=d) : centrs isa Vector{Int} ? samples[centrs] : centrs
    labels = fill(-1, length(samples))
    clusters = [Int[] for _ in centroids]
    iterations = Tuple{Vector{PhylogeneticTree}, Vector{Int64}}[]
    iteration = 1
    loss = inf
    while true
        old_labels = labels
        # Re-assign samples to clusters
        labels = [argmin(d(s, c) for c in centroids) for s in samples]
        push!(iterations, (centroids, labels))
        old_labels == labels && break
        new_loss = sum(d(s, centroids[labels[i]]) for (i, s) in enumerate(samples))
        @printf "Iteration: %2d;   ∑ₛd(s,c(s)) = %-10.2f\n" iteration new_loss
        # Compute new centroids
        empty!.(clusters)
        for (i, l) in enumerate(labels)
            push!(clusters[l], i)
        end
        println([sum(d(samples[i], centroid) for i in cluster) for (cluster, centroid) in zip(clusters, centroids)])
        centroids = [
            length(cluster) == 0 ? centroid : tropical_median_consensus(samples[cluster])
            for (cluster, centroid) in zip(clusters, centroids)
        ]
        println([sum(d(samples[i], centroid) for i in cluster) for (cluster, centroid) in zip(clusters, centroids)])
        new_loss = sum(d(s, centroids[labels[i]]) for (i, s) in enumerate(samples))
        @printf "Iteration: %2d.5; ∑ₛd(s,c(s)) = %-10.2f\n" iteration new_loss
        global L = iterations
        @assert new_loss <= loss "Loss did not decrease."
        loss = new_loss
        iteration += 1
    end
    iterations
end
iterations = cluster(samples, 5);

# Tools for plotting binary trees on four leaves
# ==============================================
""" Returns an integer that describes the topology of the tree, according to the following table:
```
    1    ((,),(,))
    2    (,(,(,)))
    3    (,((,),))
    4    ((,(,)),)
    5    (((,),),)

```
"""
function treetype(t::PhylogeneticTree)
    v = vech(t)
    if all(isapprox.([v[2], v[3], v[4], v[5]], 2.0)) # "(,),(,)"
        return 1
    elseif all(isapprox.([v[1], v[2], v[3]], 2.0))
        if isapprox(v[4], v[5]) # ",(,(,))"
            return 2
        elseif isapprox(v[5], v[6]) # ",((,),)"
            return 3
        end
    elseif all(isapprox.([v[3], v[5], v[6]], 2.0))
        if isapprox(v[1], v[2]) # "(,(,)),"
            return 4
        elseif isapprox(v[2], v[4]) # "((,),),"
            return 5
        end
    end
    nothing
end

# Assuming all trees have height one, a tree of each type can be represented by two coordinates.
# The following functions extract the coordinates, depending on type:
coordinate_functions = [
    m -> [1-m[1,2],      1-m[3,4]], # b, b''
    m -> [1-m[2,4], m[2,4]-m[3,4]], # a', b''
    m -> [1-m[2,4], m[2,4]-m[2,3]], # a', b'
    m -> [1-m[1,3], m[1,3]-m[2,3]], # a, b'
    m -> [1-m[1,3], m[1,3]-m[1,2]]  # a, b
] .∘ (m -> 1/2 * m) .∘ cophenetic_matrix

""" Returns the coordinates of `tr` in the local coordinate system of the cell corresponding to the topology type of `tr`."""
function tree_coordinates(tr::PhylogeneticTree)
    t = treetype(tr)
    return isnothing(t) ? nothing : coordinate_functions[t](tr)
end

treetype_coordinates(t) = isnothing(treetype(t)) ? nothing : (treetype(t), tree_coordinates(t))

@assert all(all(tree_coordinates(s) .>= 0) for s in samples) "Local coordinates must be nonnegative."

canvas_coordinates(t) = sum(canvas_bases[treetype(t)] .* tree_coordinates(t))

tree_from_coordinate_functions = Base.Fix1(phylogenetic_tree, Float64) .∘ [
    (b, b′) -> "(t1:$(1-b),t2:$(1-b)):$b,(t3:$(1-b′),t4:$(1-b′)):$b′;",
    (a, b) -> "t1:1.0,(t2:$(1-a),(t3:$(1-a-b),t4:$(1-a-b)):$b):$a;",
    (a, b) -> "t1:1.0,((t2:$(1-a-b),t3:$(1-a-b)):$b,t4:$(1-a)):$a;",
    (a, b) -> "(t1:$(1-a),(t2:$(1-a-b),t3:$(1-a-b)):$b):$a,t4:1.0;",
    (a, b) -> "((t1:$(1-a-b),t2:$(1-a-b)):$b,t3:$(1-a)):$a,t4:1.0;"
]
""" Convert coordinates in the local coordinate system of a cell to a PhylogeneticTree."""
tree_from_coordinates(type::Int, a, b) = tree_from_coordinate_functions[type](a, b)


const canvas_bases = [exp.(π*im.*b) for b in [
        [1//4, 3//4],   # b,b''
        [9//8, 3//4],   # a',b''
        [9//8, 6//4],   # a', b'
        [15//8, 6//4],  # a, b'
        [15//8, 1//4]   # a, b
]]

function plot_clustering(clusters, p)
    centroids, labels = clusters
    grid = [
        [0, exp(1//4*π*im)],
        [0, exp(3//4*π*im)],
        [0, exp(9//8*π*im)],
        [0, exp(6//4*π*im)],
        [0, exp(15//8*π*im)],
        [exp(1//4*π*im), exp(1//4*π*im) + exp(3//4*π*im), exp(3//4*π*im), exp(9//8*π*im), exp(6//4*π*im), exp(15//8*π*im), exp(1//4*π*im)]
    ]
    plot!(p, aspect_ratio=:equal, label="", xlims=(-1, 1), yllabelsims=(-1, 1.5))
    plot!(p, real.(grid), imag.(grid), color=:black, label="")

    s = treetype.(samples)
    c = treetype.(centroids)
    ms = (!isnothing).(s) 
    mc = (!isnothing).(c)
    s = canvas_coordinates.(samples[ms])
    c = canvas_coordinates.(centroids[mc])
    l = labels[ms]

    scatter!(p, real.(s), imag.(s), color=1 .+ l       , palette=col, marker=:circle, markersize=4, group=l)
    scatter!(p, real.(c), imag.(c), color=2:1+length(c), palette=col, marker=:star,   markersize=8, label="")
    p
end

for clusters in iterations
    display(plot_clustering(clusters, plot()))
end

xlim = (-1,1)
ylim = (-1,1.5)
res = 100

canvas_to_grid((x,y)) = (
    Int(ceil((x - xlim[1]) / (xlim[2] - xlim[1]) * res)),
    Int(ceil((y - ylim[1]) / (ylim[2] - ylim[1]) * res))
)
grid_to_canvas(i, j) = (
    xlim[1] + (i - 1) / res,
    ylim[1] + (j - 1) / res
)

inv_basis((b1, b2)) = inv([real(b1) real(b2); imag(b1) imag(b2)])

""" Convert coordinates in the canvas to coordinates in a cell, if the coordinates lie in the cell `t`.
    If the coordinates do not lie in any cell, return `nothing`."""
function canvas_to_treetype_coordinates(x, y)    
    # Find the cell in which the projection of (x,y) in local coordinates is nonnegative.
    for (t, b) in enumerate(canvas_bases)
        u, v = inv_basis(b) * [x; y]
        if ((t == 1 && 0 <= u <= 1 && 0 <= v <= 1)
          ||(t  > 1 && 0 <= u && 0 <= v && u + v <= 1))
            return (t, (u, v))
        end
    end
    return nothing
end

grid_of_trees = Array{Union{Vector{Float64}, Nothing}}(nothing, Int.(ceil.(res .* (xlim[2] - xlim[1], ylim[2] - ylim[1]))))
for i in axes(grid_of_trees, 1), j in axes(grid_of_trees, 2)
    x, y = grid_to_canvas(i, j)
    r = canvas_to_treetype_coordinates(x, y)
    isnothing(r) && continue
    (t, (u, v)) = r
    grid_of_trees[i, j] = vech(tree_from_coordinates(t, u, v))
end


xs = xlim[1] .+ (0:size(grid_of_trees, 1)-1) ./ res
ys = ylim[1] .+ (0:size(grid_of_trees, 2)-1) ./ res

function plot_clusters(clusters)
    centroids, _ = clusters
    centroids = vech.(centroids)
    hm = [
        isnothing(s) ? 0 : argmin(d(s, c) for c in centroids) 
        for s in grid_of_trees
    ]
    p = heatmap(xs, ys, transpose(hm), aspect_ratio=:equal, color=col, cmin=-.5, clim=(-.5, 9.5), alpha=.2, size=(600, 600), colorbar=false)
    hm = [
        isnothing(s) ? 0 : minimum(d(s, c) for c in centroids)
        for s in grid_of_trees
    ]
    contour!(p, xs, ys, transpose(hm), aspect_ratio=:equal, levels=50, color=col0, alpha=.5)
    plot_clustering(clusters, p)
    p
end

display.(plot_clusters.(iterations))



# %% Symthetic data equally distributed on the canvas
samples = [
    begin
        t, (a, b) = first(
            u for u in (
                canvas_to_treetype_coordinates(
                    xlim[1] + (xlim[2]-xlim[1]) * rand(),
                    ylim[1] + (ylim[2]-ylim[1]) * rand()
                )
                for _ in zip())
            if !isnothing(u)
        )
        tree_from_coordinates(t, a, b)               
    end
    for _ in 1:100
]
iterations = cluster(samples, 5);
display.(plot_clusters.(iterations));

# %% Synthetic clusters
samples = [
    begin
        t = rand(1:5)
        a, b = t == 1 ? (.5 .+ .5 .* rand(2)) : first((a, b) for (a, b) in (.3 .+ .6 .* rand(2) for _ in zip()) if a + b < 1.0)
        tree_from_coordinates(t, a, b)
    end
    for _ in 1:100
]
iterations = cluster(samples, 5);
display.(plot_clusters.(iterations));

samples = [
    begin
        t = rand([1, 3, 4])
        a, b = t == 1 ? (.5 .+ .5 .* rand(2)) : first((a, b) for (a, b) in (rand(2) for _ in zip()) if a + b < 1.0)
        tree_from_coordinates(t, a, b)
    end
    for _ in 1:100
]
iterations = cluster(samples, 2);
display.(plot_clusters.(iterations));



# %% Apicomplexa
leaq(a,b) = (a <= b) || isapprox(a, b)
is_ultrametric(m::Matrix{Float64}) = all(leaq(m[i,j], max(m[i,k], m[j,k])) for i in 1:size(m,1), j in 1:size(m,1), k in 1:size(m,1))
is_ultrametric(m::Matrix{QQFieldElem}) = all(<=(m[i,j], max(m[i,k], m[j,k])) for i in 1:size(m,1), j in 1:size(m,1), k in 1:size(m,1))

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

samples = (open("R-Data/apicomplexa.txt")
     |> readlines
    .|> (s -> phylogenetic_tree(QQFieldElem, s))
    .|> make_equidistant
)
make_equidistant(samples[1])
Random.seed!(3)
iterations = cluster(samples, 9);

# %% Consensus tree bug
f(m, t) = phylogenetic_tree(m, t)
t, s = load("consensus_tree_bug.json")
t, s = f(t...), [f(a, b) for (a, b) in s]
@assert length(unique(taxa.(s))) == 1 "All trees must have the same taxa.false"
s = samples[l]
m = tropical_median_consensus(s)
sum(d.(s, Ref(t)))
sum(d.(s, Ref(m)))

# Build and solve the corresponding LP directly, without using `tropical_median_consensus`.
V = transpose(stack(vech.(s)))
m, n = size(V)
hDiff = sum.([V[i,:]/n for i in 1:m])
for i in 1:m
    for j in 1:n
        V[i,j] -= hDiff[i]
    end
end
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
is_feasible(E)
LP = linear_program(E, l; convention=:min)
r, tx = solve_lp(LP)
t = tx[1:m]
x = tx[m+1:end]
tree = phylogenetic_tree(vech_to_matrix(x), taxa(s[1]))
sum(d.(s, Ref(tree)))

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