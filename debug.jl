using Oscar, Combinatorics

include("clustering.jl")

# %% Consensus tree bug
;f(m, t) = phylogenetic_tree(m, t)
t, s = load("consensus_tree_bug.json")
t, s = f(t...), [f(a, b) for (a, b) in s]
@assert length(unique(taxa.(s))) == 1 "All trees must have the same taxa."
tmc = tropical_median_consensus(s)
sum(d.(s, Ref(t)))
sum(d.(s, Ref(tmc)))

# Build and solve the corresponding LP directly, without using `tropical_median_consensus`.
V = transpose(stack(vech.(s[1:28])))
m, n = size(V)
hDiff = sum.([V[i,:]/n for i in 1:m])
for i in 1:m
    for j in 1:n
        V[i,j] -= hDiff[i]
    end
end
M = zeros(m*n+2, m+n)
for i in 1:m, j in 1:n
    M[(i-1) * n + j, i] = 1.0
    M[(i-1) * n + j, m + j] = 1.0
end
M[end-1, m+1:end] .=  1
M[end  , m+1:end] .= -1
v = [reshape(transpose(V), m*n); 0; 0]
E = polyhedron(-M, -v)
l = vcat(n*ones(m), zeros(n))
is_feasible(E)
LP = linear_program(E, l; convention=:min)
r, tx = solve_lp(LP)
t = tx[1:m]
sum(t)
x = tx[m+1:end]
tree = phylogenetic_tree(vech_to_matrix(x), taxa(s[1]))
sum(d.(s[1:28], Ref(tree)))

# reduce sample size
# Note that samples of cardinality divisible by 28 mostly fail. This does not seem to happen for subsets of different sizes.
for i in 1:84-27
    a, b = i, 27+i
    tmc = tropical_median_consensus(s[a:b])

    V = transpose(stack(vech.(s[a:b])))
    m, n = size(V)
    M = zeros(m*n+2, m+n)
    for i in 1:m, j in 1:n
        M[(i-1) * n + j, i] = 1.0
        M[(i-1) * n + j, m + j] = 1.0
    end
    M[end-1, m+1:end] .=  1
    M[end  , m+1:end] .= -1
    v = [reshape(transpose(V), m*n); 0; 0]
    E = polyhedron(-M, -v)
    l = vcat(n*ones(m), zeros(n))
    LP = linear_program(E, l; convention=:min)
    r, tx = solve_lp(LP)
    t = tx[1:m]
    x = tx[m+1:end]
    tree = phylogenetic_tree(vech_to_matrix(x), taxa(s[1]))


    if sum(d.(s[a:b], Ref(tmc))) - sum(d.(s[a:b], Ref(tree))) > 1
        println(i)
    end
end

# reconstruct polymake side for a,b = 1,28
# Note that if we compare the resulting cophenetic matrices
cophenetic_matrix(tree) 
cophenetic_matrix(tmc)
# particularly the entries [6,8] and [7,8] in the c. matrix of tmc are significantly bigger than those of tree (even after subtracting the constant factor that Andrei adds)

S = zeros(28,28) # essentially V 
for k in 1:m
    for (i,j) in combinations(1:8,2)
            S[k,8*(i-1) - Int(i*(i+1)//2)+ j] = cophenetic_matrix(s[k])[i,j]
    end
end
S
# Map row vectors of S to hyperplane H 
hDiff = sum.([S[i,:]/n for i in 1:m])
for i in 1:m
    for j in 1:n
        S[i,j] -= hDiff[i]
    end
end
S 

# As expected, when comparing the result of the tropical_median function with x above the last two entries diverge significantly after making up for the 
tm = Vector(Polymake.call_function(:tropical,:tropical_median, S))
filter(i -> abs(tm[i] - x[i]) > 1, 1:28)

# debugging tropical_median
# 1. computing tropical vertices
m, n = 28, 28
supply = n*ones(m); demand = m*ones(n)
flowMatrix = Polymake.call_function(:graph,:optimal_transport_plan, -S, supply, demand)

# Check correctness of flowMatrix
trans_poly = transportation_polytope(supply, demand)
convert(Vector{Rational}, vcat(v[1,:],v[2,:]))
trans_LP = linear_program(trans_poly, matrix(QQFieldElem,V))

trop_vert = Matrix(Polymake.call_function(:tropical, :facets_matrix, S, flowMatrix))
# 2. computing the average of the tropical vertices
tm = zeros(n)
r, c = size(trop_vert)
for i in 1:c
    tm[i] = 1/r*(sum(trop_vert[:,i]))
end
# The result is similar to the result of the tropical_median function in polymake (up to adding some multiple of the all ones vector)
# The problem must be in the computation of the tropical vertices then

fw_set = Polymake.call_function(:tropical, :fw_set, S)
dim(polyhedron(fw_set))

function facets_matrix(costMatrix, flowMatrix)
    m = nrows(costMatrix)
    n = ncols(costMatrix)
   
    @assert m == nrows(flowMatrix) && n == ncols(flowMatrix) "The dimension of the cost matrix does not match the dimension of the flow matrix."
   
    g = Graph{Undirected}(m+n)
   
    # here we prefer the coordinate nodes to have the indices from 0 to $n - 1
    for i in 1:m
        for j in 1:n
            # if ($flowMatrix -> elem($i, $j) != 0) {
            if flowMatrix[i,j] != 0
                add_edge!(g, j, i + n);
            end
        end
    end
   
    nv = n_vertices(g)
    multip = zeros(nv) # dual prices
    labels = (-1)*ones(Int, nv) # represent the connecting components

    d = 0

    for j in 1:n
        if labels[j] == -1
            d += 1
            dfs_initialize_mult(g, j, d, labels, multip, costMatrix, n);
        end
    end
    println(labels)
	
    delta = zeros(d, d)
    bigM = 0
    for i in 1:m 
        for j in 1:n
            if bigM < multip[n+i] + multip[j] - costMatrix[i, j] 
                bigM = multip[n+i] + multip[j] - costMatrix[i, j]
            end
        end
    end

    for i in 1:d
        for j in 1:d 
            if (i != j)
                delta[i,j] = bigM;
            end
        end
    end
   
    for i in 1:m, j in 1:n
        diff = multip[n + i] + multip[j] - costMatrix[i, j]
        if delta[labels[n + i], labels[j]] > diff
            delta[labels[n + i], labels[j]] = diff
        end
    end
	
	# Floyd-Warshall algorithm
    for k in 1:d, i in 1:d, j in 1:d
        if delta[i, j] > delta[i, k] + delta[k, j]
            delta[i, j] = delta[i, k] + delta[k, j]
        end
    end

	# facets -> (i,j) contains the maximal value of x_i - x_j over the polytrope
    facets = [multip[i] - multip[j] + delta[labels[i], labels[j]] for i in 1:n, j in 1:n]

    return facets
end

function dfs_initialize_mult(G, node, lab, labList, mult, edgeCosts, N)
    labList[node] = lab
    # println("G: ", G)
    # println("node: ", node)
    # println("N: ", N)

    for next in neighbors(G, node)
        if labList[next] == -1
            i, j = 1,1
            if node < N + 1
                j = node
                i = next - N
            else
                i = node - N
                j = next
            end

            mult[next] = edgeCosts[i, j] - mult[node]
            dfs_initialize_mult(G, next, lab, labList, mult, edgeCosts, N)
        end
    end
end

flowMatrix = Polymake.call_function(:graph,:optimal_transport_plan, -S, supply, demand)
facets = facets_matrix(S, Matrix(flowMatrix))
n_edges(g)