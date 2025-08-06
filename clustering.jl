using Oscar
import Oscar: PhylogeneticTree

using LinearAlgebra
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
function vech(x::PhylogeneticTree)
    vech(cophenetic_matrix(x))
end

function d(x::PhylogeneticTree, y::PhylogeneticTree)
    nx = length(taxa(x))
    ny = length(taxa(y))
    @assert nx == ny "Cannot compare trees with different numbers of leaves."
    vx = vech(x)
    vy = vech(y)
    sum(vy .- vx) - nx*minimum(vy .- vx)
end

using JSON
join(readlines(open("/homes/combi/lenzen/ClusteringPhyloTrees/R-Data/newick-100-4.json")))
JSON.parse(ans)
samples = phylogenetic_tree.(Float64, map(v -> v[1], ans))

iteration = 1
centroids = samples[1:3] 
labels = [-1 for _ in samples]
clusters = rand(samples, 3)
while true
    println("Iteration: ", iteration)
    old_labels = labels
    labels = [argmin(d(s, t) for t in centroids) for s in samples]
    if(old_labels == labels)
        break
    end
    clusters = [Int[] for _ in centroids]
    for (i, l) in enumerate(labels)
        push!(clusters[l], i)
    end
    for t in eachindex(centroids)
        length(clusters[t]) == 0 && continue
        centroids[t] = tropical_median_consensus(samples[clusters[t]])
    end
    iteration += 1
end

using Plots
bases_args = [
    [1//4, 3//4],
    [3//4, 9//8],
    [9//8, 6//4],
    [6//4, 15//8],
    [15//8, 1//4]
]
types = [
    "((,),(,))",
    "(,(,(,)))",
    "(,((,),))",
    "((,(,)),)",
    "(((,),),)",
]
coordinates = [
    m -> (1-m[1,2], 1-m[3,4]),
    m -> (1-m[2,3], m[2,3]-m[3,4]),
    
]
re = r"[t0-9.:;]+"
gettype(t::PhylogeneticTree) = replace(newick(t), re => "")

bases = [exp.(π*im.*b) for b in bases_args]
p = plot(aspect_ratio=:equal, legend=:none, xlims=(-1.5, 1.5), ylims=(-1.5, 1.5))
for (b, _) in bases
    plot!(p, [0, real(b)], [0, imag(b)], label="", color=:black)
end

types = findfirst.(isequal.(gettype.(samples)), Ref(types))
coordinates = 
