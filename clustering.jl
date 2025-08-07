using Oscar, LinearAlgebra, Plots, JSON, OrderedCollections, Random
import Oscar: PhylogeneticTree

""" Convert a square matrix to the vector of its upper triangular part, in column-major order."""
function vech(A::AbstractMatrix{T}) where T
    m = LinearAlgebra.checksquare(A)
    v = Vector{T}(undef, (m*(m+1))>>1)
    k = 0
    for j = 1:m, i = j:m
        @inbounds v[k += 1] = A[i,j]
    end
    return v
end

""" Convert a PhylogeneticTree to a vector in ℝᵉ/𝟏ℝ, where e = choose(#leaves, 2)."""
vech(t::PhylogeneticTree) = vech(cophenetic_matrix(t))

function d(x::PhylogeneticTree, y::PhylogeneticTree)
    nx = length(taxa(x))
    ny = length(taxa(y))
    @assert nx == ny "Cannot compare trees with different numbers of leaves."
    vx = vech(x)
    vy = vech(y)
    sum(vy .- vx) - nx*minimum(vy .- vx)
end

head(x) = x[1]

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
    .|> head
    .|> strip_parentheses
    .|> (do_normalize ? normalize : identity)
    .|> Base.Fix1(phylogenetic_tree, Float64)
)

samples = read_newick_json("R-Data/newick-100-4.json", true)

Random.seed!(2)
function cluster(samples, k)
    iteration = 1
    centroids = rand(samples, k)
    labels = [-1 for _ in samples]
    clusters = [Int[] for _ in centroids]
    iterations = Tuple{Vector{PhylogeneticTree}, Vector{Int64}}[]
    while true
        old_labels = labels
        labels = [argmin(d(s, t) for t in centroids) for s in samples]
        if(old_labels == labels)
            break
        end
        println("Iteration: ", iteration)
        for c in clusters
            empty!(c)
        end
        for (i, l) in enumerate(labels)
            push!(clusters[l], i)
        end
        for t in eachindex(centroids)
            length(clusters[t]) == 0 && continue
            centroids[t] = tropical_median_consensus(samples[clusters[t]])
        end
        push!(iterations, (centroids, labels))
        iteration += 1
    end
    iterations
end
iterations = cluster(samples, 5)

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
treetype(t::PhylogeneticTree) = findfirst(isequal(break_ties(typestring(t))), types)

samples_per_type = mergewith(vcat, Dict{Int64, Vector{Int64}}(), [Dict(treetype(t) => [i]) for (i, t) in enumerate(samples)]...)
# Assuming all trees have height one, a tree of each type can be represented by two coordinates.
# The following functions extract the coordinates, depending on type:
coordinate_functions = [
    m -> [1-m[1,2],      1-m[3,4]], # b, b''
    m -> [1-m[2,4], m[2,4]-m[3,4]], # a', b''
    m -> [1-m[2,4], m[2,4]-m[2,3]], # a', b'
    m -> [1-m[1,3], m[1,3]-m[2,3]], # a, b'
    m -> [1-m[1,3], m[1,3]-m[1,2]]  # a, b
] .∘ (t -> .5*cophenetic_matrix(t))

tree_from_coordinates = phylogenetic_tree .∘ [
    (b, b′′) -> 2 .* [0 1-b 1 1; 1-b 0 1 1; 1 1 0 1-b′′; 1 1 1-b′′ 0],
    (a′, b′′) -> 2 .* [0 1 1 1; 1 0 1-a′ 1-a′; 1 1-a′ 0 1-a′-b′′; 1 1-a′ 1-a′-b′′ 0],
    # (a′, b′) -> 2 .* [0 1 1 1; 
]

coordinates(t) = coordinate_functions[treetype(t)](t)
treetypes = treetype.(samples)
local_coords = hcat(coordinates.(samples)...)
@assert all(local_coords .>= 0) "Coordinates must be non-negative."

function plot_clustering(clusters)
    centroids, labels = clusters
    bases = [exp.(π*im.*b) for b in [
            [1//4, 3//4],   # b,b''
            [9//8, 3//4],   # a',b''
            [9//8, 6//4],   # a', b'
            [15//8, 6//4],  # a, b'
            [15//8, 1//4]   # a, b
    ]]

    grid = [
        [0, exp(1//4*π*im)],
        [0, exp(3//4*π*im)],
        [0, exp(9//8*π*im)],
        [0, exp(6//4*π*im)],
        [0, exp(15//8*π*im)],
        [exp(1//4*π*im), exp(1//4*π*im) + exp(3//4*π*im), exp(3//4*π*im), exp(9//8*π*im), exp(6//4*π*im), exp(15//8*π*im), exp(1//4*π*im)]
    ]
    p = plot(aspect_ratio=:equal, label="", xlims=(-1, 1), ylims=(-1, 1.5))
    for b in grid
        plot!(p, real.(grid), imag.(grid), color=:black, label="")
    end
    canvas_coordinates = [sum(coordinates(s) .* bases[treetype(s)]) for s in samples]
    scatter!(p, real.(canvas_coordinates), imag.(canvas_coordinates), color=labels)
    centroids_coordinates = [sum(coordinates(s) .* bases[treetype(s)]) for s in centroids]
    scatter!(p, real.(centroids_coordinates), imag.(centroids_coordinates), color=1:length(centroids_coordinates), marker=:star, markersize=10, label="Centroids")
end

p = nothing
for clusters in iterations
    display(plot_clustering(clusters))
end

xlim = (-1,1)
ylim = (-1,1.5)
res = 100
grid = Vector{Int64}(undef, res .* Int.(ceil.((xlim[2] - xlim[1]) .* (ylim[2] - ylim[1]))))
for i, j in eachindex(grid)
    x = xlim[1] + (i - 1) * (xlim[2] - xlim[1]) / res
    y = ylim[1] + (j - 1) * (ylim[2] - ylim[1]) / res
    grid[i, j] = argmin(d(PhylogeneticTree(x, y), t) for t in centroids)
end