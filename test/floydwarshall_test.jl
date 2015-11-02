function floydwarshall_test()
    A = load_matrix_network("all_shortest_paths_example")
    A.nzval += abs(minimum(A)) + 1 # remove negative edges
    m = size(A,1)
    D2 = zeros(Float64,m,m)
    P2 = zeros(Float64,m,m)
    for i=1:m
        (d,p) = dijkstra(A,i)
        D2[i,:] = d
        P2[i,:] = p
    end
    (D,P) = floydwarshall(A)
    
    if !isequal(D,D2)
        error("Floyd Warshall failed: Incorrect distances")
    end
    
    if !isequal(P,P2)
        error("Floyd Warshall failed: Incorrect predecessors")
    end
    return true
end