n = 2
a = 5
sites = transpose(hcat([
    [[0,0,0]];
    [
        [n+a/n, -a/n, -n]
        for n in 1:15
    ]
]...))
sites = sites[:, 2:end] .- sites[:,1]

sites = [
    1 1;
    -1 -1
]

N = 1001
f = 100
offset = -751

d(a,b) = sum(b .- a) - 3*min(0, minimum(b .- a))
# d(a,b) = sum(b .- a) + 3*max(0, maximum(a .- b))

grid = [
    argmin(d(((i,j) .+ offset) ./ f, s) for s in eachrow(sites))
    for i in 1:N, j in 1:N
]

using Plots
r = ((1:N) .+ offset) ./ f
default(colormap=:viridis)
heatmap(r, r, grid, aspect_ratio=:equal)
scatter!(eachcol(sites)..., zcolor=1:size(sites, 1), markerstrokecolor=:white)