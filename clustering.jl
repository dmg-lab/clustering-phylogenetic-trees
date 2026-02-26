using Oscar, LinearAlgebra, Plots, JSON, OrderedCollections, Random, StatsBase, Printf, Combinatorics
import Oscar: PhylogeneticTree

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

vech(A::QQMatrix) = vech(collect(A))

""" Convert a vector to a symmetric matrix, in row-major order, with diagonals set to zero."""
vech_to_matrix(v::AbstractVector) = vech_to_matrix_internal(v)

vech_to_matrix(v::AbstractVector{QQFieldElem}) = matrix(vech_to_matrix_internal(v))

function vech_to_matrix_internal(v::AbstractVector{T}) where T
    n = Int((sqrt(1 + 8*length(v)) + 1) ÷ 2)
    @assert n * (n - 1) ÷ 2 == length(v) "Vector length must be n(n-1)/2 for some n."
    A = zeros(T, n, n)
    k = 0
    for i in 1:n, j in i+1:n
        @inbounds A[j, i] = A[i, j] = v[k += 1]
    end
    return A
end


""" Convert a PhylogeneticTree to a vector in ℝᵉ/𝟏ℝ, where e = choose(#leaves, 2)."""
vech(t::PhylogeneticTree) = vech(cophenetic_matrix(t))

""" Asymmetric tropical distance. """
function d_trop(x::Vector, y::Vector)
    sum(y .- x) - length(x)*minimum(y .- x)
end

function d_trop(x::PhylogeneticTree, y::PhylogeneticTree)
    nx = length(taxa(x))
    ny = length(taxa(y))
    @assert nx == ny "Cannot compare trees with different numbers of leaves."
    d_trop(vech(x), vech(y))
end

function normalize(newick_string)
    i = 0
    r(_) = "t$(i += 1)"
    replace(newick_string, r"t[0-9]+" => r)
end

strip_parentheses(newick_string) = replace(newick_string, r"^\((.*)\);$" => s"\1;")

""" Reads `filename` as a json file containing a list of newick strings. """
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

"""
    farthest_point_sampling_rand(samples, k; d=d_trop)

Chooses k centroids according to the (randomized) k-means++ initialization scheme.
"""
function farthest_point_sampling_rand(samples, k; d=d_trop)
    cs = Vector{PhylogeneticTree}(undef, k)
    cs[1] = rand(samples)
    ds = d.(Ref(cs[1]), samples)
    for i in 2:k
        cs[i] = sample(samples, pweights(ds))
        ds = min.(ds, d.(Ref(cs[i]), samples))
    end
    return cs
end

"""
    farthest_point_sampling_strict(samples, k; d=d_trop)

Chooses a random point from `samples`, and chooses the `k` points that are farthest away
from all other previously selected points.
"""
function farthest_point_sampling_strict(samples, k; d=d_trop)
    cs = Vector{PhylogeneticTree}(undef, k)
    cs[1] = rand(samples)
    ds = d.(Ref(cs[1]), samples)
    for i in 2:k
        cs[i] = samples[argmax(ds)]
        ds = min.(ds, d.(Ref(cs[i]), samples))
    end
    return cs
end

loss(centroids, labels, samples; d=d_trop) = sum(d(s, centroids[labels[i]]) for (i, s) in enumerate(samples))

""" 
    cluster(samples::Vector{T}, k::Union{Int, Vector{Int}, Vector{T}}; d=d_trop, median_func=tropical_median_consensus) where T

k-means-clustering w.r.t. the distance function `d`. Initialization of the centroids is as follows: If `k` is
* an `Int`: select initial centroids by k-means++-clustering (see `farthest_point_sampling_rand()`)
* a `Vector{Int}`: use the samples indexed by `k` as initial centroids
* a `Vector` of the same type as `samples`: use these centroids.

Returns:

A vector representing the state after each iteration of the algorithm
The vector consists of pairs `(centroids, labels)`, where `c::Vector{T}` is the list of centroids after the respective iteration,
and `labels::Vector{Int}` is a vector such that `samples[i]` belongs to the `lables[i]`-th cluster.
"""
function cluster(samples::Vector{T}, k::Union{Int, Vector{Int}, Vector{T}}; d=d_trop, median_func=tropical_median_consensus) where T
    # farthest-point sampling of initial centroids
    centroids = k isa Int ? farthest_point_sampling_rand(samples, k; d=d) : k isa Vector{Int} ? samples[k] : k
    
    labels = fill(-1, length(samples))      # index i of the cluster each sample s belongs to
    clusters = [Int[] for _ in centroids]   # samples s belonging to cluster i
    iterations = Tuple{Vector{T}, Vector{Int64}}[]
    while true
        # Re-assign samples to clusters
        old_labels = labels
        labels = [argmin(d(s, c) for c in centroids) for s in samples]
        push!(iterations, (centroids, labels))

        # break if clustering is stationary
        old_labels == labels && break

        # Compute new centroids
        empty!.(clusters)
        for (i, l) in enumerate(labels)
            push!(clusters[l], i)
        end
        centroids = [
            length(cluster) == 0 ? centroid : median_func(samples[cluster])
            for (cluster, centroid) in zip(clusters, centroids)
        ]
    end
    return iterations
end
