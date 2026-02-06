using Oscar

function aws(V, w, ws)
    Ew = polyhedron(-V,w)
    Ews = polyhedron(-V, ws)
    result = findmin(x -> 
        findmax(y -> 
            findmin(i -> 
                (dot(V[i,:],Vector(x)) + w[i]) / (dot(V[i,:],Vector(y)) + ws[i]), filter(i -> dot(V[i,:],Vector(y)) + ws[i] != 0, 1:n_rows(V)))[1]
        ,vertices(Ews))[1],
    vertices(Ew))[1]
    return result
end

function non_split_prime_part(n,w)
    h = hypersimplex(2,n)
    V = matrix(vertices(h))
    h = h.pm_polytope
    s = Polymake.polytope.splits(h.VERTICES,h.GRAPH.ADJACENCY,h.FACETS,n-1)
    sub = subdivision_of_points(V, w)
    s_in_sub = Polymake.polytope.splits_in_subdivision(h.VERTICES, sub.pm_subdivision.MAXIMAL_CELLS, s)
    
    w_minus_w0 = zeros(nrows(V))
    s = matrix(QQ, s)[:,2:end]
    for j in s_in_sub
        row = s[j+1,:]
        ws = V*row
        for i in 1:length(ws)
            if ws[i] < 0 ws[i] = 0 end
        end
        a = aws(V, w, ws)
        w_minus_w0 += a * ws
    end
    return w_minus_w0
end 

n = 5;
w = Vector([1,2,3,4,3,4,2,1,3//2,3//2]);
w_minus_w0 = non_split_prime_part(n, w);
V = matrix(vertices(hypersimplex(2,n)))
maximal_cells(subdivision_of_points(V, w_minus_w0))