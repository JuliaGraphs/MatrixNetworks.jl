function dijkstra_test()
    (A,xy,labels) = load_matrix_network_metadata("airports")
    A = -A; # fix funny encoding of airport data
    lax = 247; rst = 355
    (d,pred) = dijkstra(A,lax)
    return true
end