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
function cluster(samples, centrs::Union{Int, Vector{Int}, Vector{PhylogeneticTree}}; d=d, median_func=tropical_median_consensus)
    # farthest-point sampling of initial centroids
    centroids = centrs isa Int ? farthest_point_sampling_rand(samples, centrs; d=d) : centrs isa Vector{Int} ? samples[centrs] : centrs
    
    labels = fill(-1, length(samples))      # index i of the cluster each sample s belongs to
    clusters = [Int[] for _ in centroids]   # samples s belonging to cluster i
    iterations = Tuple{Vector{PhylogeneticTree}, Vector{Int64}}[]
    loss = inf                              # loss with current clustering
    while true
        # Re-assign samples to clusters
        old_labels = labels
        labels = [argmin(d(s, c) for c in centroids) for s in samples]
        push!(iterations, (centroids, labels))
        # @assert loss <= (loss = sum(d(s, centroids[labels[i]]) for (i, s) in enumerate(samples))) "Loss did not decrease."

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
        # @assert loss <= (loss = sum(d(s, centroids[labels[i]]) for (i, s) in enumerate(samples))) "Loss did not decrease."
    end
    iterations
end
