include("clustering.jl")
include("plotting.jl")

# %% Read some phylogenetic trees and cluster into k=5 clusters.
samples = read_newick_json("example-data/newick-100-4.json", true)
iterations = cluster(samples, 5);

# %% Plot the clustering iterations
for (i,f) in enumerate(iterations)
    p = plot_clusters(samples, f)
    plot!(p, legend=false, size=(250,300), axis=false, bottom_margin=-10mm, top_margin=-10mm, left_margin=-10mm, right_margin=-10mm)
    display(p)
    # savefig(p, "clustering_iteration_$(i).pdf")
end

# Generate uniformly distributed random trees on four leaves and cluster
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
for i in iterations
    display(plot_clusters(samples, i))
end

#%% Generate the data uniformly taken from the δ-coarse topologies.
#    We will find that the clusters usually align with the clusters one "sees".
Random.seed!(10) # for this seed (and all other parameters unchanged), there will be four iterations until the clustering converges.
δ = 0.15
samples = PhylogeneticTree{Float64}[]
while length(samples) <= 200
    t = rand(1:5)
    a, b = rand(2)
    if (t == 1 && a >= δ && b >= δ) || (t != 1 && a >= δ && a+b <= 1)
        push!(samples, tree_from_coordinates(t, a, b))
    end
end
display(plot_trees(samples))
# Cluster, and plot (and save) the clustering
iterations = cluster(samples, 3)
for (i,f) in enumerate(iterations)
    p = plot_clusters(samples, f)
    display(p)
    # savefig(p, "clustering_iteration_$(i).pdf")
end

# %% Generate synthetic data clustered by the five topologies
#    We will find that the clusters usually don't align with the clusters one sees.
Random.seed!(5)
δ = 0.2
samples = PhylogeneticTree{Float64}[]
while length(samples) <= 200
    t = rand(1:5)
    a, b = rand(2)
    if a >= δ && b >= δ && (t == 1 || a+b <= 1)
        push!(samples, tree_from_coordinates(t, a, b))
    end
end

iterations = cluster(samples, 3);
for (i,f) in enumerate(iterations)
    p = plot_clusters(samples, f)
    display(p)
    # savefig(p, "clustering_iteration_$(i).pdf")
end

iterations = cluster(samples, 5);
for (i,f) in enumerate(iterations)
    p = plot_clusters(samples, f)
    display(p)
    # savefig(p, "clustering_iteration_$(i).pdf")
end