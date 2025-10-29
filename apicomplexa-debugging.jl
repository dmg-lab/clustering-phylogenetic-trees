# More debugging stuff (Lena)
# ===========================
# reduce sample size
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
for k in 1:8
    for (i,j) in combinations(1:8,2)
            S[k,8*(i-1) - Int(i*(i+1)//2)+ j] = cophenetic_matrix(s[k])[i,j]
    end
end
# As expected, when comparing the result of the tropical_median function with x above the last two entries diverge significantly after making up for the
tm = Vector(Polymake.call_function(:tropical,:tropical_median, S))
filter(i -> abs(tm[i] - x[i]) > 1, 1:28)

# debugging tropical_median
# 1. computing tropical vertices
m, n = 28, 28
supply = n*ones(m); demand = m*ones(n)
flowMatrix = Polymake.call_function(:graph,:optimal_transport_plan, -S, supply, demand)
trop_vert = Polymake.call_function(:tropical, :facets_matrix, S, flowMatrix)
# 2. computing the average of the tropical vertices
tm = zeros(n)
r, c = size(trop_vert)
for i in 1:c
    tm[i] = 1/r*(sum(trop_vert[:,i]))
end
# The result is similar to the result of the tropical_median function in polymake (up to adding some multiple of the all ones vector)
# The problem must be in the computation of the tropical vertices then