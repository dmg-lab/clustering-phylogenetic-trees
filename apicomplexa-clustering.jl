using Combinatorics, ProgressMeter
include("apicomplexa-functions.jl")

# Read data and make all trees equidistant by extending edges to leaves.
samples = (open("../data/R-data/apicomplexa.txt")
     |> readlines
    .|> (s -> phylogenetic_tree(Float64, s)) # <== also with Float64 and QQFieldElem
    .|> make_equidistant
)

#%% Run k-clusterings with different seeds for the (randomized) centroid initialization.
#   We want to find the seed which minimizes the total loss.
#   CAUTION! The following loop will take quite a while.
#   See the code cell below for the seed that minimizes the loss.
k = 17
losses = @showprogress map(1:100) do s
    Random.seed!(s)
    try
        centroids, labels = cluster(samples, k; median_func=tropical_median_consensus_div)[end]
        loss(centroids, labels, samples)
    catch e
        @warn "cluster failed for seed $s: $e"
        Inf
    end   
end
loss_min, seed_min = findmin(losses)
h_losses = histogram(losses, bins = range(loss_min-1, ceil(maximum(losses))+10, step=20), legend = false, color = :gray, linecolor = :black, grid = true, xticks = floor(minimum(losses)):20:ceil(maximum(losses))+5, xtickfontsize = 14, ytickfontsize = 14) # savefig("losses.pdf")

#%% We ran this before. The following seed attains the minimum.
Random.seed!(31)
centroids, labels = cluster(samples, k; median_func=tropical_median_consensus_div)[end]
loss(centroids, labels, samples)

# To make sense of the centroids as trees, make all edge lengths positive and all trees have same height.
# Since d_trop is invariant under denormalization, the loss stays the same:
centroids = denormalize.(centroids)
centroids = (centroids .+ (maximum(height.(centroids)) .- height.(centroids)))
loss(centroids, labels, samples)
Plots.plot(plot_phylo.(centroids)...)
# using Plots.PlotMeasures
# for (i,q) in enumerate(centroids)
#     p = plot_phylo(q)
#     plot!(p, size=(150, 150), bottom_margin=-6mm, top_margin=-1mm, left_margin=-5mm, right_margin=-1mm)
#     display(p)
#     # savefig(p, "$i.pdf")
# end

#%% Understand better the nature of the different clusters (See paper)
samples_per_cluster = let
    r = [Int[] for _ in centroids]
    for (i, l) in enumerate(labels)
        push!(r[l], i)
    end
    r
end
# The clusters have various sizes.
length.(samples_per_cluster)

len_SPC = length.(samples_per_cluster) 
lonely = findall(i -> i == 1, len_SPC)
lonelies_15 = Dict{Int, Vector{Int}}([i => [] for i in lonely]...)
lonelies_17 = Dict{Int, Vector{Int}}([i => [] for i in lonely]...)
lonelies_10 = Dict{Int, Vector{Int}}([i => [] for i in lonely]...)

k = 10
@showprogress for s in 1:100
    Random.seed!(s)
    try
        cen, lab = cluster(samples, k; median_func=tropical_median_consensus_div)[end]
        loss(cen, lab, samples)
        
        spc = let
            r = [Int[] for _ in cen]
            for (i, l) in enumerate(lab)
                push!(r[l], i)
            end
            r
        end
 
        lonelies = findall(i -> i == 1, length.(spc))

        for i in lonelies
            j = findfirst(t -> isapprox(d(cen[i], centroids_0[t]), 0), lonely)
            if !isnothing(j)
                push!(lonelies_10[lonely[j]], s)
            end
        end

    catch e
        @warn "cluster failed for seed $s: $e"
        Inf
    end   
end

# diameter per cluster
diameter = zeros(k)
for i in 1:k
    s = samples[samples_per_cluster[i]]
    if length(s) != 1
        pairs = collect(Combinatorics.combinations(1:length(s), 2))
        distances = [d(s[i], s[j]) for (i,j) in pairs]
        diameter[i] = maximum(distances)
    end
end

# "density"
density = [diameter[i]/length(samples_per_cluster[i]) for i in 1:k]

# ambiguity in timing of split (Pf, Pv)
ambi_PfPv = [(maximum(cophenetic_matrix(centroids[i])) - cophenetic_matrix(centroids[i])[4,5])/2 for i in 1:k]
ambi_BbTa = [(maximum(cophenetic_matrix(centroids[i])) - cophenetic_matrix(centroids[i])[1,6])/2 for i in 1:k]
ambi_EtTg = [(maximum(cophenetic_matrix(centroids[i])) - cophenetic_matrix(centroids[i])[3,7])/2 for i in 1:k]

histogram(ambi_PfPv, bins=range(0,maximum(ambi_PfPv)+0.1, step=0.125), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14, ylims=(0,5))
savefig("pictures/ambiguity_PfPv.pdf")
histogram(ambi_BbTa, bins=range(0,maximum(ambi_PfPv)+0.1, step=0.125), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_BbTa.pdf")
histogram(ambi_EtTg, bins=range(0,maximum(ambi_EtTg)+0.1, step=0.125), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_EtTg.pdf")

for i in 1:17
    if !isempty(filter(M -> all(isapprox.(M,cophenetic_matrix(centroids[i]))), cophenetic_matrix.(samples)))
        println(i)
    end
end

# counts (Bb, Ta) (-> 1rst and 6th taxon) in samples
# counts (Pf, Pv) (-> 4th and 5th taxon) in samples
count(cophenetic_matrix.(samples)) do M
    row = M[4,:]
    isapprox(row[5], minimum(row[[1:3;5:end]])) && !isapprox(row[5], maximum(row))
end

# counts (Et, Tg) in centroids
c = 0
for ci in 1:17
    row = cophenetic_matrix(centroids[ci])[3,:]
    if isapprox(row[7], minimum(row[[1:2;4:end]])) && !isapprox(row[7], maximum(row))
        c += len_SPC[ci]
    end
end

histogram([(maximum(M[1,:]) - M[1,6])/2 for M in cophenetic_matrix.(samples)], bins=range(0,8,step=.25), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14, ylims=(0,150))
savefig("pictures/ambiguity_samples_BbTa.pdf")
histogram([(maximum(M[3,:]) - M[3,7])/2 for M in cophenetic_matrix.(samples)], bins=range(0,8,step=.25), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14, ylims=(0,150))
savefig("pictures/ambiguity_samples_EtTg.pdf")
histogram([(maximum(M[4,:]) - M[4,5])/2 for M in cophenetic_matrix.(samples)], bins=range(0,8,step=.25), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14,  ylims=(0,150))
savefig("pictures/ambiguity_samples_PfPv.pdf")