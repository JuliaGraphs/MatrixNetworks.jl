function bipartite_matching_test()
    W = sprand(10,8,0.5)
    bipartite_matching(W)
    bipartite_matching([10;12;13],[1;2;3],[3;2;4])
    return true
end