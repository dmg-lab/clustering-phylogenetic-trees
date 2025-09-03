include("clustering.jl")
# Tools for plotting binary trees on four leaves
# ==============================================
""" Returns an integer that describes the topology of the tree, according to the following table:
```
    1    ((,),(,))
    2    (,(,(,)))
    3    (,((,),))
    4    ((,(,)),)
    5    (((,),),)

```
"""
function treetype(t::PhylogeneticTree)
    v = vech(t)
    if all(isapprox.([v[2], v[3], v[4], v[5]], 2.0)) # "(,),(,)"
        return 1
    elseif all(isapprox.([v[1], v[2], v[3]], 2.0))
        if isapprox(v[4], v[5]) # ",(,(,))"
            return 2
        elseif isapprox(v[5], v[6]) # ",((,),)"
            return 3
        end
    elseif all(isapprox.([v[3], v[5], v[6]], 2.0))
        if isapprox(v[1], v[2]) # "(,(,)),"
            return 4
        elseif isapprox(v[2], v[4]) # "((,),),"
            return 5
        end
    end
    nothing
end

# Assuming all trees have height one, a tree of each type can be represented by two coordinates.
# The following functions extract the coordinates, depending on type:
coordinate_functions = [
    m -> [1-m[1,2],      1-m[3,4]], # b, b''
    m -> [1-m[2,4], m[2,4]-m[3,4]], # a', b''
    m -> [1-m[2,4], m[2,4]-m[2,3]], # a', b'
    m -> [1-m[1,3], m[1,3]-m[2,3]], # a, b'
    m -> [1-m[1,3], m[1,3]-m[1,2]]  # a, b
] .∘ (m -> 1/2 * m) .∘ cophenetic_matrix

""" Returns the coordinates of `tr` in the local coordinate system of the cell corresponding to the topology type of `tr`."""
function tree_coordinates(tr::PhylogeneticTree)
    t = treetype(tr)
    return isnothing(t) ? nothing : coordinate_functions[t](tr)
end

treetype_coordinates(t) = isnothing(treetype(t)) ? nothing : (treetype(t), tree_coordinates(t))

@assert all(all(tree_coordinates(s) .>= 0) for s in samples) "Local coordinates must be nonnegative."

canvas_coordinates(t) = sum(canvas_bases[treetype(t)] .* tree_coordinates(t))

tree_from_coordinate_functions = Base.Fix1(phylogenetic_tree, Float64) .∘ [
    (b, b′) -> "(t1:$(1-b),t2:$(1-b)):$b,(t3:$(1-b′),t4:$(1-b′)):$b′;",
    (a, b) -> "t1:1.0,(t2:$(1-a),(t3:$(1-a-b),t4:$(1-a-b)):$b):$a;",
    (a, b) -> "t1:1.0,((t2:$(1-a-b),t3:$(1-a-b)):$b,t4:$(1-a)):$a;",
    (a, b) -> "(t1:$(1-a),(t2:$(1-a-b),t3:$(1-a-b)):$b):$a,t4:1.0;",
    (a, b) -> "((t1:$(1-a-b),t2:$(1-a-b)):$b,t3:$(1-a)):$a,t4:1.0;"
]
""" Convert coordinates in the local coordinate system of a cell to a PhylogeneticTree."""
tree_from_coordinates(type::Int, a, b) = tree_from_coordinate_functions[type](a, b)


const canvas_bases = [exp.(π*im.*b) for b in [
        [1//4, 3//4],   # b,b''
        [9//8, 3//4],   # a',b''
        [9//8, 6//4],   # a', b'
        [15//8, 6//4],  # a, b'
        [15//8, 1//4]   # a, b
]]

function plot_clustering(clusters, p)
    centroids, labels = clusters
    grid = [
        [0, exp(1//4*π*im)],
        [0, exp(3//4*π*im)],
        [0, exp(9//8*π*im)],
        [0, exp(6//4*π*im)],
        [0, exp(15//8*π*im)],
        [exp(1//4*π*im), exp(1//4*π*im) + exp(3//4*π*im), exp(3//4*π*im), exp(9//8*π*im), exp(6//4*π*im), exp(15//8*π*im), exp(1//4*π*im)]
    ]
    plot!(p, aspect_ratio=:equal, label="", xlims=(-1, 1), yllabelsims=(-1, 1.5))
    plot!(p, real.(grid), imag.(grid), color=:black, label="")

    s = treetype.(samples)
    c = treetype.(centroids)
    ms = (!isnothing).(s) 
    mc = (!isnothing).(c)
    s = canvas_coordinates.(samples[ms])
    c = canvas_coordinates.(centroids[mc])
    l = labels[ms]

    scatter!(p, real.(s), imag.(s), color=1 .+ l       , palette=col, marker=:circle, markersize=4, group=l)
    scatter!(p, real.(c), imag.(c), color=2:1+length(c), palette=col, marker=:star,   markersize=8, label="")
    p
end

for clusters in iterations
    display(plot_clustering(clusters, plot()))
end

xlim = (-1,1)
ylim = (-1,1.5)
res = 100

canvas_to_grid((x,y)) = (
    Int(ceil((x - xlim[1]) / (xlim[2] - xlim[1]) * res)),
    Int(ceil((y - ylim[1]) / (ylim[2] - ylim[1]) * res))
)
grid_to_canvas(i, j) = (
    xlim[1] + (i - 1) / res,
    ylim[1] + (j - 1) / res
)

inv_basis((b1, b2)) = inv([real(b1) real(b2); imag(b1) imag(b2)])

""" Convert coordinates in the canvas to coordinates in a cell, if the coordinates lie in the cell `t`.
    If the coordinates do not lie in any cell, return `nothing`."""
function canvas_to_treetype_coordinates(x, y)    
    # Find the cell in which the projection of (x,y) in local coordinates is nonnegative.
    for (t, b) in enumerate(canvas_bases)
        u, v = inv_basis(b) * [x; y]
        if ((t == 1 && 0 <= u <= 1 && 0 <= v <= 1)
          ||(t  > 1 && 0 <= u && 0 <= v && u + v <= 1))
            return (t, (u, v))
        end
    end
    return nothing
end

grid_of_trees = Array{Union{Vector{Float64}, Nothing}}(nothing, Int.(ceil.(res .* (xlim[2] - xlim[1], ylim[2] - ylim[1]))))

for i in axes(grid_of_trees, 1), j in axes(grid_of_trees, 2)
    x, y = grid_to_canvas(i, j)
    r = canvas_to_treetype_coordinates(x, y)
    isnothing(r) && continue
    (t, (u, v)) = r
    grid_of_trees[i, j] = vech(tree_from_coordinates(t, u, v))
end

xs = xlim[1] .+ (0:size(grid_of_trees, 1)-1) ./ res
ys = ylim[1] .+ (0:size(grid_of_trees, 2)-1) ./ res

function plot_clusters(clusters)
    centroids, _ = clusters
    centroids = vech.(centroids)
    hm = [
        isnothing(s) ? 0 : argmin(d(s, c) for c in centroids) 
        for s in grid_of_trees
    ]
    p = heatmap(xs, ys, transpose(hm), aspect_ratio=:equal, color=col, cmin=-.5, clim=(-.5, 9.5), alpha=.2, size=(600, 600), colorbar=false)
    hm = [
        isnothing(s) ? 0 : minimum(d(s, c) for c in centroids)
        for s in grid_of_trees
    ]
    contour!(p, xs, ys, transpose(hm), aspect_ratio=:equal, levels=50, color=col0, alpha=.5)
    plot_clustering(clusters, p)
    p
end

