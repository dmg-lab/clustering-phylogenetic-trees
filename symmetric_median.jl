###############################################################################
#  Symmetric tropical Fermat–Weber set via min-cost flow  (Float64 + Oscar)
#
#  Data: rows v_1..v_m of V are points in R^n / R·1. The FW set is
#     argmin_x  Σ_i d_tr(x, v_i),   d_tr(x, v) = max_j (x_j - v_j) - min_j (x_j - v_j).
#
#  Idea: n · Σ_i d_tr(x, v_i) is the LP dual of a min-cost flow on the
#  network I -> J -> K (I, K: one node per data point, J: one per coordinate):
#     arcs  i -> j  with cost  v_ij,   every i supplies n,
#     arcs  j -> k  with cost -v_kj,   every k demands  n.
#  The optimal dual potentials on J are exactly the FW points. They are the
#  potentials π with π_b - π_a ≤ cost(a,b) on every arc of the residual graph
#  of an optimal flow, so the FW set is the polytrope
#     { x : x_k - x_j ≤ K[j,k] },
#  where K is the J-block of the residual graph's all-pairs shortest paths
#  (its tropical Kleene star). Equivalently it is the min-tropical convex
#  hull of the rows of K.
###############################################################################

using Oscar

"""
    tfw(V; tol=nothing) -> (K, Fmin)

Rows of `V` are the data points. Returns the n×n Kleene star `K` of the
symmetric tropical Fermat–Weber set  { x : x_k - x_j ≤ K[j,k] }  (the
min-tropical convex hull of the rows of `K`), and `Fmin = min Σ_i d_tr(x, v_i)`.
`tol` (default 1e-9·max(1, max|V|)) only absorbs rounding.
"""
function tfw(V; tol=nothing)
    V = Matrix{Float64}(V)
    all(isfinite, V) || throw(ArgumentError("±Inf/NaN entries are not supported"))
    m, n = size(V)
    tol = something(tol, 1e-9 * max(1.0, maximum(abs, V)))

    # Network: nodes I = 1:m, J = m+1:m+n, K = m+n+1:2m+n.
    arcs = vcat([(i, m + j, V[i, j]) for i in 1:m for j in 1:n],          # i -> j
                [(m + j, m + n + k, -V[k, j]) for k in 1:m for j in 1:n])  # j -> k
    demand = [fill(-1.0n, m); zeros(n); fill(1.0n, m)]   # inflow - outflow
    N, E = 2m + n, length(arcs)

    # Min-cost flow as an LP over Float64 with Oscar. The balance equations
    # sum to 0, so one (node 1) is dropped to keep them consistent in floats.
    A = zeros(N, E)
    for (e, (u, v, _)) in enumerate(arcs)
        A[u, e] -= 1
        A[v, e] += 1
    end
    nonneg = (-Float64.((1:E) .== (1:E)'), zeros(E))     # f ≥ 0
    P = polyhedron(Float64, nonneg, (A[2:N, :], demand[2:N]))
    cost, f = solve_lp(linear_program(P, last.(arcs); convention=:min))

    # Residual graph of the optimal flow (arcs are uncapacitated): every arc,
    # plus the reversed arc with negated cost for every arc carrying flow.
    # Floyd–Warshall gives its tropical Kleene star.
    D = fill(Inf, N, N)
    for v in 1:N
        D[v, v] = 0.0
    end
    for ((u, v, c), fe) in zip(arcs, f)
        D[u, v] = min(D[u, v], c)
        fe > tol && (D[v, u] = min(D[v, u], -c))
    end
    for r in 1:N, a in 1:N, b in 1:N
        D[a, b] = min(D[a, b], D[a, r] + D[r, b])
    end
    minimum(D[v, v] for v in 1:N) >= -tol * N || error("negative cycle in residual graph")

    K = D[m+1:m+n, m+1:m+n]
    K .-= [K[j, j] for j in 1:n]       # diag(K) = 0 (absorbs rounding)
    return K, -Float64(cost) / n
end

"Tropical distance d_tr(x, v) = max_j (x_j - v_j) - min_j (x_j - v_j)."
dtr(x, v) = maximum(x .- v) - minimum(x .- v)

"Σ_i d_tr(x, v_i): the FW objective; equals `Fmin` exactly on the FW set."
fw_objective(V, x) = sum(dtr(x, V[i, :]) for i in axes(V, 1))

"Whether `x` lies in the FW set described by `K`."
in_fw_set(K, x; tol=1e-9) = all(x' .- x .<= K .+ tol)    # [j,k] = x_k - x_j

"""
    fw_polytope(K; tol=1e-9)

The FW set as an Oscar tropical polyhedron (min convention): the tropical
convex hull of the rows of `K`. Oscar's tropical polyhedra are rational, so
each entry is replaced by the simplest rational within `tol`.
"""
fw_polytope(K; tol=1e-9) =
    tropical_convex_hull(matrix(QQ, [QQ(rationalize(BigInt, x; tol=tol)) for x in K]), min)

# ----------------------------------------------------------------------------
# Example
# ----------------------------------------------------------------------------
V = [0.0 0.0 0.0; 0.0 3.0 1.0; 0.0 1.0 4.0; 0.0 2.0 2.0]
K, Fmin = tfw(V)                                   # Fmin = 7
[fw_objective(V, K[j, :]) for j in 1:3] .≈ Fmin    # all true
P = fw_polytope(K); vertices(P)

function sym_tropical_median(V)

end