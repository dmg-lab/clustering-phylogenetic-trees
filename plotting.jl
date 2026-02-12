include("clustering.jl")
const col0 = palette([:black])
const col1 = palette([colorant"#ffffffff", palette(:tab10)...])

"""
    treetype(t::PhylogeneticTree)

Returns an integer that describes the topology of the tree, according to the following table:
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

"""
    tree_coordinates(tr::PhylogeneticTree)

Returns the coordinates of `tr` in the local coordinate system of the cell corresponding to the topology type of `tr`.
"""
function tree_coordinates(tr::PhylogeneticTree)
    coordinate_functions = [
        m -> [1-m[1,2],      1-m[3,4]], # b, b''
        m -> [1-m[2,4], m[2,4]-m[3,4]], # a', b''
        m -> [1-m[2,4], m[2,4]-m[2,3]], # a', b'
        m -> [1-m[1,3], m[1,3]-m[2,3]], # a, b'
        m -> [1-m[1,3], m[1,3]-m[1,2]]  # a, b
    ] .∘ (m -> 1/2 * m) .∘ cophenetic_matrix
    t = treetype(tr)
    return isnothing(t) ? nothing : coordinate_functions[t](tr)
end

treetype_coordinates(t) = isnothing(treetype(t)) ? nothing : (treetype(t), tree_coordinates(t))

canvas_coordinates(t) = sum(canvas_bases[treetype(t)] .* tree_coordinates(t))

tree_from_coordinate_functions = Base.Fix1(phylogenetic_tree, Float64) .∘ 
""" Convert coordinates in the local coordinate system of a cell to a PhylogeneticTree."""
function tree_from_coordinates(type::Int, a, b)
    @assert (0. <= a <= 1. && 0. <= b <= 1. && (type == 1 || a+b <= 1.)) "Coordinates do not specify tree of given type."
    phylogenetic_tree(Float64, [
        (b, b′) -> "(t1:$(1-b),t2:$(1-b)):$b,(t3:$(1-b′),t4:$(1-b′)):$b′;",
        (a, b) -> "t1:1.0,(t2:$(1-a),(t3:$(1-a-b),t4:$(1-a-b)):$b):$a;",
        (a, b) -> "t1:1.0,((t2:$(1-a-b),t3:$(1-a-b)):$b,t4:$(1-a)):$a;",
        (a, b) -> "(t1:$(1-a),(t2:$(1-a-b),t3:$(1-a-b)):$b):$a,t4:1.0;",
        (a, b) -> "((t1:$(1-a-b),t2:$(1-a-b)):$b,t3:$(1-a)):$a,t4:1.0;"
    ][type](a,b))
end

const canvas_bases = [exp.(π*im.*b) for b in [
        [1//4, 3//4],   # b,b''
        [9//8, 3//4],   # a',b''
        [9//8, 6//4],   # a', b'
        [15//8, 6//4],  # a, b'
        [15//8, 1//4]   # a, b
]]

function plot_trees!(p, samples)
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
    ms = (!isnothing).(s) 
    s = canvas_coordinates.(samples[ms])

    scatter!(p, real.(s), imag.(s),  color=:black, marker=:circle, markersize=4, legend=false)
    p
end
plot_trees(samples, args...) = plot_trees!(plot(legend=false, size=(250,300), aspect_ratio=:equal, axis=false, bottom_margin=-10mm, top_margin=-10mm, left_margin=-10mm, right_margin=-10mm, args...), samples)

function plot_clustering!(samples, clusters, p)
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

    scatter!(p, real.(s), imag.(s), color=0 .+ l     , palette=col1[2:end], marker=:circle, markersize=4, group=l)
    scatter!(p, real.(c), imag.(c), color=1:length(c), palette=col1[2:end], marker=:star,   markersize=8, label="")
    p
end


const xlim = (-1,1)
const ylim = (-1,1.5)
const res = 100

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

const grid_of_trees = let
    grid_of_trees = Array{Union{Vector{Float64}, Nothing}}(nothing, Int.(ceil.(res .* (xlim[2] - xlim[1], ylim[2] - ylim[1]))))
    for i in axes(grid_of_trees, 1), j in axes(grid_of_trees, 2)
        x, y = grid_to_canvas(i, j)
        r = canvas_to_treetype_coordinates(x, y)
        isnothing(r) && continue
        (t, (u, v)) = r
        grid_of_trees[i, j] = vech(invokelatest(tree_from_coordinates, t, u, v))
    end
    grid_of_trees
end

const xs = xlim[1] .+ (0:size(grid_of_trees, 1)-1) ./ res
const ys = ylim[1] .+ (0:size(grid_of_trees, 2)-1) ./ res

function plot_clusters!(p, samples, clusters)
    centroids, _ = clusters
    centroids = vech.(centroids)
    hm = [
        isnothing(s) ? 0 : argmin(d(s, c) for c in centroids) 
        for s in grid_of_trees
    ]
    p = heatmap!(xs, ys, transpose(hm), aspect_ratio=:equal, color=col1, cmin=-.5, clim=(-.5, 9.5), alpha=.2)
    hm = [
        isnothing(s) ? 0 : minimum(d(s, c) for c in centroids)
        for s in grid_of_trees
    ]
    contour!(p, xs, ys, transpose(hm), aspect_ratio=:equal, levels=50, color=col0, alpha=.5)
    plot_clustering!(samples, clusters, p)
    p
end

plot_clusters(samples, clusters, args...) = plot_clusters!(
    plot(legend=false, size=(250,300), aspect_ratio=:equal, axis=false,
         bottom_margin=-10mm, top_margin=-10mm, left_margin=-10mm, right_margin=-10mm, args...),
    samples, clusters
)

function explain_visualization(from=(1,0,0))
    grid = [
        exp(1//4*π*im),
        exp(3//4*π*im),
        exp(9//8*π*im),
        exp(6//4*π*im),
        exp(15//8*π*im),
    ]
    O = vech(tree_from_coordinates(from...))
    p = plot(aspect_ratio=:equal, size=(800,800), showaxis=false, legend=false, xlims=(-1.0,1.0), ylims=(-1.0,1.5))
    for (i, g) in enumerate(grid)
        plot!(p, [0, real(g)], [0, imag.(g)], color=:black, label="", aspect_ratio=:equal, size=(800,800))
        # annotate!(p, 1.1*real(g), 1.1 * imag(g), "$i")
    end
    g = [exp(1//4*π*im), exp(1//4*π*im) + exp(3//4*π*im), exp(3//4*π*im), exp(9//8*π*im), exp(6//4*π*im), exp(15//8*π*im), exp(1//4*π*im)]
    plot!(p, real.(g), imag.(g), color=:black)
    hm1 = [
        isnothing(s) ? 0 : d(O, s)
        for s in grid_of_trees
    ]
    hm2 = [
        isnothing(s) ? 0 : d(s, O)
        for s in grid_of_trees
    ]
    lvls = Int(ceil(max(maximum(hm1), maximum(hm2))))
    contour!(p, xs, ys, transpose(hm1), aspect_ratio=:equal, levels=lvls, clims=(0,lvls+1), color=:black)
    contour!(p, xs, ys, transpose(hm2), aspect_ratio=:equal, levels=lvls, clims=(0,lvls+1), color=:red, colorbar=false)
    return p
end
# explain_visualization((1,.0, .0))
# savefig(explain_visualization(), "visualization_explanation.pdf")
