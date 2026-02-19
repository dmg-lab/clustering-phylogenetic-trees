using Combinatorics

include("apicomplexa-functions.jl")

samples = (open("data/R-data/apicomplexa.txt")
     |> readlines
    #  |> Base.Fix2(getindex, load("R-Data/apicomplexa_path_subset.txt")) # Subset that causes the problem to appear
    .|> (s -> phylogenetic_tree(Float64, s)) # <== also with Float64 and QQFieldElem
    .|> make_equidistant
    # .|> normalize_cophenetic_matrix
)

"""
    tropical_median_consensus_div(samples)

If `length(samples)` is divisible by the dimension,
then (for numeric instability reasons), one sample is excluded,
before `invoking tropical_median_consensus2`.
"""
function tropical_median_consensus_div(samples)
    n = length(taxa(samples[1]))
    if length(samples) % (n*(n-1)//2) == 0
        println("⚠️ Excluding one sample for divisibility reasons.")
        r = Int(rand(UInt) % length(samples)) + 1
        tropical_median_consensus2(samples[[1:r; r+2:end]])
    else
        tropical_median_consensus2(samples[1:end])
    end
end

"""
    repeat_until_success(f, args...; kwargs...)

Runs `f(args...; kwargs...)` until no exception is raised.
Of course, this only makes sense for non-deterministic functions.
"""
function repeat_until_success(f, args...; kwargs...)
    while true
        try
            return f(args...; kwargs...)
        catch e
            println("⚠️ Failed; try again")
            continue
        end
    end
end

using ProgressMeter
ProgressMeter.ncalls(::typeof(argmin), ::Function, arg) = length(arg)
ProgressMeter.ncalls(::typeof(findmin), ::Function, arg) = length(arg)

# The following loop will take a while, which is why it is commented out.
# The seed that minimizes the loss of the clustering 
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
h_losses = histogram(losses, 
    bins = range(loss_min-1, ceil(maximum(losses))+10, step=20),
    legend = false,
    color = :gray,
    linecolor = :black,
    grid = true,
    xticks = floor(minimum(losses)):20:ceil(maximum(losses))+5,
    xtickfontsize = 14,
    ytickfontsize = 14)
savefig("losses.pdf")

Random.seed!(31)
# Random.seed!(s[2])
centroids, labels = cluster(samples, k; median_func=tropical_median_consensus_div)[end]
loss(centroids, labels, samples)

# Make centroids have all positive edge lengths and all of same height
centroids = denormalize.(centroids)
centroids = (centroids .+ (maximum(height.(centroids)) .- height.(centroids)))
loss(centroids, labels, samples)

Plots.plot(plot_phylo.(centroids)...)
using Plots.PlotMeasures
for (i,q) in enumerate(centroids)
    p = plot_phylo(q)
    plot!(p, size=(150, 150), bottom_margin=-6mm, top_margin=-1mm, left_margin=-5mm, right_margin=-1mm)
    display(p)
    # savefig(p, "$i.pdf")
end

samples_per_cluster = let
    r = [Int[] for _ in centroids]
    for (i, l) in enumerate(labels)
        push!(r[l], i)
    end
    r
end
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

histogram(ambi_PfPv, bins=range(0,maximum(ambi_PfPv)+0.1, step=0.001), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_PfPv.pdf")
histogram(ambi_BbTa, bins=range(0,maximum(ambi_BbTa)+0.1, step=0.001), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_BbTa.pdf")
histogram(ambi_EtTg, bins=range(0,maximum(ambi_EtTg)+0.1, step=0.001), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_EtTg.pdf")

for i in 1:17
    if !isempty(filter(M -> all(isapprox.(M,cophenetic_matrix(centroids[i]))), cophenetic_matrix.(samples)))
        println(i)
    end
end

# counts (Et, Tg) in samples
count(cophenetic_matrix.(samples)) do M
    row = M[6,:]
    isapprox(row[1], minimum(row[[1:5;7:end]])) && !isapprox(row[1], maximum(row))
end

# counts (Et, Tg) in centroids
c = 0
for ci in 1:17
    row = cophenetic_matrix(centroids[ci])[3,:]
    if isapprox(row[7], minimum(row[[1:2;4:end]])) && !isapprox(row[7], maximum(row))
        c += len_SPC[ci]
    end
end

histogram([(maximum(M[1,:]) - M[1,6])/2 for M in cophenetic_matrix.(samples)], bins=range(0,10,step=.3), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_samples_BbTa.pdf")
histogram([(maximum(M[3,:]) - M[3,7])/2 for M in cophenetic_matrix.(samples)], bins=range(0,10,step=.3), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_samples_EtTg.pdf")
histogram([(maximum(M[4,:]) - M[4,5])/2 for M in cophenetic_matrix.(samples)], bins=range(0,10,step=.3), legend=false, color=:gray, xtickfontsize = 14,ytickfontsize = 14)
savefig("ambiguity_samples_PfPv.pdf")