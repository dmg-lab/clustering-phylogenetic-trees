include("apicomplexa-functions.jl")

samples = (open("R-Data/apicomplexa.txt")
     |> readlines
    #  |> Base.Fix2(getindex, load("R-Data/apicomplexa_path_subset.txt")) # Subset that causes the problem to appear
    .|> (s -> phylogenetic_tree(QQFieldElem, s)) # <== also with Float64 and QQFieldElem
    .|> make_equidistant
    .|> normalize_cophenetic_matrix
)

tm1 = tropical_median_consensus(samples)
tm2 = tropical_median_consensus2(samples)
println(d(tm1, tm2))
println(d(tm2, tm1))