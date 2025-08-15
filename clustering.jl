using Oscar, LinearAlgebra, Plots, JSON, OrderedCollections, Random, StatsBase 
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

Random.seed!(2)

function find_centroids(samples, k; d=d)
    cs = Vector{PhylogeneticTree}(undef, k)
    cs[1] = rand(samples)
    ds = d.(Ref(cs[1]), samples)
    for i in 2:k
        cs[i] = sample(samples, pweights(ds))
        ds = min.(ds, d.(Ref(cs[i]), samples))
    end
    return cs
end

""" k-means-clustering w.r.t. the distance function `d`.
    Either provide `centroids` as a vector of indices into `samples`, or an integer.
    In the latter case, the function samples `centroids` many centroids at random."""
function cluster(samples, centroids::Union{Int, Vector{Int}}; d=d)
    cs = centroids isa Int ? find_centroids(samples, centroids; d=d) : samples[centroids]
    labels = fill(-1, length(samples))
    clusters = [Int[] for _ in cs]
    iterations = Tuple{Vector{PhylogeneticTree}, Vector{Int64}}[]
    iteration = 1
    while true
        old_labels = labels
        labels = [argmin(d(s, c) for c in cs) for s in samples]
        push!(iterations, (cs, labels))
        if old_labels == labels
            break
        end
        print("Iteration: ", iteration)
        for c in clusters
            empty!(c)
        end
        for (i, l) in enumerate(labels)
            push!(clusters[l], i)
        end
        for t in eachindex(cs)
            length(clusters[t]) == 0 && continue
            cs[t] = tropical_median_consensus(samples[clusters[t]])
        end
        println("; d = ", sum(d(s, cs[labels[i]]) for (i, s) in enumerate(samples)))
        iteration += 1
    end
    iterations
end
iterations = cluster(samples, 5);

# Tools for plotting binary trees on four leaves
# ==============================================

# There are only the following five binary trees on four leaves:
types = [
    "(,),(,)",
    ",(,(,))",
    ",((,),)",
    "(,(,)),",
    "((,),),",
]

# Extract the type string from a PhylogeneticTree, as in `types`:
typestring(t::PhylogeneticTree) = replace(newick(t), r"[t0-9.:;e-]+" => "")

# Break ties in non-binary newick string by adding parentheses around commas:
function break_ties(newick_string)
    d = [false]
    e = [""]
    o = ""
    for c in newick_string
        if c == '('
            o *= c
            push!(d, false)
            push!(e, "")
        elseif c == ')'
            o *= c * pop!(e)
            pop!(d)
        elseif c == ','
            if !d[end]
                o *= c
                d[end] = true
            else
                o *= "(,"
                push!(d, true)
                e[end] *= ")"
            end
        end
    end
    o *= pop!(e)
    @assert isempty(e) "Unmatched parentheses"
    o
end
# Extract the index of the type:
# treetype(t::PhylogeneticTree) = findfirst(isequal(break_ties(typestring                                    
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

samples_per_type = mergewith(vcat, Dict{Int64, Vector{Int64}}(), [Dict(treetype(t) => [i]) for (i, t) in enumerate(samples)]...)
# Assuming all trees have height one, a tree of each type can be represented by two coordinates.
# The following functions extract the coordinates, depending on type:
coordinate_functions = [
    m -> [1-m[1,2],      1-m[3,4]], # b, b''
    m -> [1-m[2,4], m[2,4]-m[3,4]], # a', b''
    m -> [1-m[2,4], m[2,4]-m[2,3]], # a', b'
    m -> [1-m[1,3], m[1,3]-m[2,3]], # a, b'
    m -> [1-m[1,3], m[1,3]-m[1,2]]  # a, b
] .∘ (m -> 1/2 * m) .∘ cophenetic_matrix

function tree_coordinates(tr)
    t = treetype(tr)
    return isnothing(t) ? nothing : coordinate_functions[t](tr)
end

treetype_coordinates(t) = isnothing(treetype(t)) ? nothing : (treetype(t), tree_coordinates(t))

@assert all(all(tree_coordinates(s) .>= 0) for s in samples) "Local coordinates must be nonnegative."

canvas_coordinates(t) = sum(canvas_bases[treetype(t)] .* tree_coordinates(t))

tree_from_coordinates = Base.Fix1(phylogenetic_tree, Float64) .∘ [
    (b, b′) -> "(t1:$(1-b),t2:$(1-b)):$b,(t3:$(1-b′),t4:$(1-b′)):$b′;",
    (a, b) -> "t1:1.0,(t2:$(1-a),(t3:$(1-a-b),t4:$(1-a-b)):$b):$a;",
    (a, b) -> "t1:1.0,((t2:$(1-a-b),t3:$(1-a-b)):$b,t4:$(1-a)):$a;",
    (a, b) -> "(t1:$(1-a),(t2:$(1-a-b),t3:$(1-a-b)):$b):$a,t4:1.0;",
    (a, b) -> "((t1:$(1-a-b),t2:$(1-a-b)):$b,t3:$(1-a)):$a,t4:1.0;"
]

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
    plot!(p, aspect_ratio=:equal, label="", xlims=(-1, 1), ylims=(-1, 1.5))
    plot!(p, real.(grid), imag.(grid), color=:black, label="")

    s = treetype.(samples)
    c = treetype.(centroids)
    ms = (!isnothing).(s) 
    mc = (!isnothing).(c)
    s = canvas_coordinates.(samples[ms])
    c = canvas_coordinates.(centroids[mc])
    l = labels[ms]

    scatter!(p, real.(s), imag.(s), color=1 .+ l       , palette=col, marker=:circle, markersize=4, label="")
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
    grid_of_trees[i, j] = vech(tree_from_coordinates[t](u, v))
end



xs = xlim[1] .+ (0:size(grid_of_trees, 1)-1) ./ res
ys = ylim[1] .+ (0:size(grid_of_trees, 2)-1) ./ res
hm = nothing
plots = []
for clusters in iterations
    centroids, labels = clusters
    centroids = vech.(centroids)
    hm = [
        isnothing(s) ? 0 : argmin(d(s, c) for c in centroids) 
        for s in grid_of_trees
    ]
    p = heatmap(xs, ys, transpose(hm), aspect_ratio=:equal, color=col, cmin=-.5, clim=(-.5, 9.5), alpha=.2, legend=false)
    hm = [
        isnothing(s) ? 0 : minimum(d(s, c) for c in centroids)
        for s in grid_of_trees
    ]
    contour!(p, xs, ys, transpose(hm), aspect_ratio=:equal, levels=50, color=col0, alpha=.5)
    plot_clustering(clusters, p)
end
plot(plots[1:4]..., layout=(1,4), legend=false, size=(800*length(plots),800))

# %% Symthetic data
samples = PhylogeneticTree{Float64}[]
for s in eachrow(rand(100,2) .* [xlim[2]-xlim[1] ylim[2]-ylim[1]] .+ [xlim[1] ylim[1]])
       u = canvas_to_treetype_coordinates(s...)
       isnothing(u) && continue
       (t, (a, b)) = u
       push!(samples, tree_from_coordinates[t](a, b))
end


# %% Apicomplexa

""" Make a phylogenetic tree equidistant by adding sufficient lengths to the edges
    adjacent to the leaves. All lengths of interiour edges remain the same. """
function make_equidistant(tree::PhylogeneticTree)
    graph = adjacency_tree(tree)
    edge_lengths = tree.pm_ptree.EDGE_LENGTHS;
    height = tree.pm_ptree.NODE_HEIGHTS[1]
    leaf_indices = Dict(tree.pm_ptree.LEAVES[t]+1 => i for (i, t) in enumerate(taxa(tree)))
    m = cophenetic_matrix(tree)
    function f(v, h)
        for w in outneighbors(graph, v)
            h′ = h + edge_lengths[Edge(v, w)]
            if outdegree(graph, w) == 0
                j = leaf_indices[w]
                m[:,j] .+= height - h′
                m[j,:] .+= height - h′
                m[j,j]  -= 2*(height - h′)
            else
                f(w, h′)
            end
        end
    end
    f(1, 0.0)
    new_tree = phylogenetic_tree(m, taxa(tree))
    @assert is_equidistant(new_tree) "The new tree is not equidistant."
    return new_tree
end

samples = (open("R-Data/apicomplexa.txt")
     |> readlines
    .|> (s -> phylogenetic_tree(Float64, s))
    .|> make_equidistant
)
cluster(samples, 5)
