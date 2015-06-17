W = sprand(10,8,0.5)
bipartite_matching(W)
M_out = bipartite_matching([10;12;13],[1;2;3],[3;2;4])
M_out.weight
M_out.cardinality
M_out.match
MatrixNetworks.create_sparse(bipartite_matching(W)) # get the sparse matrix
MatrixNetworks.edgelist(bipartite_matching(W)) # get the edgelist

using MAT
function bipartite_matching_test()
    W = sprand(10,8,0.5)
    bipartite_matching(W)
    M_out = bipartite_matching([10;12;13],[1;2;3],[3;2;4])
    return M_out;
#     M_out.weight
#     M_out.cardinality
#     M_out.match
#     MatrixNetworks.create_sparse(bipartite_matching(W)) # get the sparse matrix
#     MatrixNetworks.edgelist(bipartite_matching(W)) # get the edgelist
end