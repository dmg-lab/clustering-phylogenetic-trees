using Combinatorics
include("apicomplexa-functions.jl")

# Read data and make all trees equidistant by extending edges to leaves.
samples = (open("data/lungfish_filtered.txt")
     |> readlines
    .|> (s -> phylogenetic_tree(Float64, s))
    .|> make_equidistant
)

#%% Run tropical k-means clustering.
k = 2
Random.seed!(1)
centroids, labels = cluster(samples, k; median_func=tropical_median_consensus_div)[end]
total_loss = loss(centroids, labels, samples)
println("k = $k, total loss = $total_loss")

samples_per_cluster = let
    r = [Int[] for _ in centroids]
    for (i, l) in enumerate(labels)
        push!(r[l], i)
    end
    r
end
println("cluster sizes: ", length.(samples_per_cluster))

# Make centroids interpretable as trees: positive edge lengths, common height.
# d_trop is invariant under denormalization, so the loss is unchanged.
centroids_readable = denormalize.(centroids)
centroids_readable = centroids_readable .+ (maximum(height.(centroids_readable)) .- height.(centroids_readable))
for (i, c) in enumerate(centroids_readable)
    println("centroid $i: ", newick(c))
end

Plots.plot(plot_phylo.(centroids)...)
