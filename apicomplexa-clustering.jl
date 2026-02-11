include("apicomplexa-functions.jl")
using Phylo

samples = (open("data/R-data/apicomplexa.txt")
     |> readlines
    #  |> Base.Fix2(getindex, load("R-Data/apicomplexa_path_subset.txt")) # Subset that causes the problem to appear
    .|> (s -> phylogenetic_tree(Float64, s)) # <== also with Float64 and QQFieldElem
    .|> make_equidistant
    # .|> normalize_cophenetic_matrix
)

function tmc(samples)
    n = length(taxa(samples[1]))
    if length(samples) % (n*(n-1)//2) == 0
        println("⚠️ Excluding one sample for divisibility reasons.")
        r = Int(rand(UInt) % length(samples)) + 1
        tropical_median_consensus2(samples[[1:r; r+2:end]])
    else
        tropical_median_consensus2(samples[1:end])
    end
end

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

ks = 2:10
ls = Tuple{Float64, Float64}[]
loss((centroids, labels)) = sum(d(s, centroids[labels[i]]) for (i, s) in enumerate(samples))
using ProgressMeter
@showprogress for k in ks
    ts = Float64[]
    for _ in 1:25
        push!(ts, loss(cluster(samples, k; median_func=tmc)[end]))
    end
    push!(ls, mean_and_std(ts))
end
using Plots
plot(ks, ls, xlabel="k", ylabel="Loss")


# p = repeat_until_success(cluster, samples, k; median_func=tmc)
# Plots.plot([Plots.plot(parsenewick("($(newick(denormalize(q))[1:end-1]));")) for q in p[end][1]]...)
# length.(filter(i -> i == j,p[end][2]) for j in 1:k)

# reduce sample size for testing
#samples = samples[1:12] # <== also try without this line
# samples = samples[1:28]

#@assert all(is_equidistant.(samples))
#@assert all(is_normalized.(samples))

# Tropical median consensus tree via Andrei's implementation
# @time mt1 = tropical_median_consensus(samples)
# @assert is_equidistant(mt1)
# is_normalized(mt1)
# nm1 = normalize_cophenetic_matrix(cophenetic_matrix(mt1))
# d1 = Float64(sum(d.(samples, Ref(mt1))))

# # Tropical median consensus tree calling tropical_median directly:
# @time mt2 = tropical_median_consensus2(samples)
# @assert is_equidistant(mt2)
# is_normalized(mt2)
# # nm2 = normalize_cophenetic_matrix(cophenetic_matrix(mt2))
# d2 = Float64(sum(d.(samples, Ref(mt2))))
# d(mt1, mt2)

# # The following is the straight forward way by solving the primal LP; while that should also be correct, it might take long time.
# @time mt3 = tropical_median_consensus3(samples)
# # @assert is_equidistant(mt3)
# # nm3 = normalize_cophenetic_matrix(cophenetic_matrix(mt3))
# # d3 = sum(d.(samples, Ref(mt3)))
# d(mt2, mt3)

