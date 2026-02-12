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


centroids, labels = cluster(samples, 17; median_func=tropical_median_consensus_div)[end]

using ProgressMeter
ProgressMeter.ncalls(::typeof(argmin), ::Function, arg) = length(arg)
centroids, labels = @showprogress argmin(((c,l),)->(loss(c,l,samples)), (cluster(samples, 17; median_func=tropical_median_consensus_div)[end] for _ in 1:100))

Plots.plot(plot_phylo.(centroids)...)
for (i,q) in enumerate(centroids)
    savefig(plot_phylo(centroids), "$i.pdf")
end