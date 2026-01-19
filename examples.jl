include("clustering.jl")
include("plotting.jl")
using Plots.PlotMeasures

# %% Read and cluster some phylogenetic trees
samples = read_newick_json("data/R-data/newick-100-4.json", true)
iterations = cluster(samples, 5);

# %% Plot the clustering iterations
display.(plot_clusters.(iterations))
for (i,f) in enumerate(iterations)
    p = plot_clusters(f)
    plot!(p, legend=false, size=(250,300), axis=false, bottom_margin=-10mm, top_margin=-10mm, left_margin=-10mm, right_margin=-10mm)
    display(p)
    savefig(p, "clustering_iteration_$(i).pdf")
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
display.(plot_clusters.(iterations));

# %% Generate synthetic data clustered around three of the five topologies
samples = [
    begin
        t = rand([1,3,4])
        a, b = t == 1 ? (.5 .+ .5 .* rand(2)) : first((a, b) for (a, b) in (.3 .+ .6 .* rand(2) for _ in zip()) if a + b < 1.0)
        tree_from_coordinates(t, a, b)
    end
    for _ in 1:100
]
iterations = cluster(samples, 5);
display.(plot_clusters.(iterations));
