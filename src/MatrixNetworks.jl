module MatrixNetworks

using Compat

if VERSION < v"0.4-"
    using Docile
end

"""
Module ``MatrixNetworks``: Documentation on the module

- Option 1: start with a sparse matrix A:
- example: ``M = MatrixNetwork(A)``
- Option 2: start with row and column indexes for the nonzeros in the matrix
- example: ``M = MatrixNetwork(ei,ej)``

Available functions: (use ?function_name to get more documentation)
- bfs
- bipartite_matching
- clustercoeffs
- corenums
- cosineknn
- csr_to_sparse
- dfs
- dijkstra
- dirclustercoeffs
- floydwarshall
- largest_component
- mst_prim
- scomponents
- sparse_to_csr

You can check the readme file here: \n
"https://github.com/nassarhuda/MatrixNetworks.jl/blob/master/README.md"
"""
MatrixNetworks

include("MatrixNetwork.jl")
export MatrixNetwork, sparse_transpose

include("scomponents.jl")
include("csr_to_sparse.jl")
include("sparse_to_csr.jl")
include("bipartite_matching.jl")
include("bfs.jl")
include("dfs.jl")
include("clustercoeffs.jl")
include("corenums.jl")
include("floydwarshall.jl")

include("manage_data.jl")
export load_matrix_network, load_matrix_network_metadata, load_matrix_network_all, 
    matrix_network_datasets
include("largest_component.jl")
include("cosineknn.jl")
include("dirclustercoeffs.jl")
include("dijkstra.jl")
include("mst_prim.jl")

include("spectral.jl")
export fiedler_vector, sweepcut, spectral_cut, bestset, SweepcutProfile

# export everything to make them accessible as functions
export bipartite_matching, edge_list, create_sparse,
bipartite_matching_setup, bipartite_matching_indicator, bfs,
dfs, clustercoeffs, corenums, scomponents, strong_components_map,
enrich, csr_to_sparse, floydwarshall, largest_component,
sparse_to_csr, cosineknn, dirclustercoeffs, dijkstra, mst_prim, mst_prim_matrix,
csr_to_sparse_matrix, edge_indicator

include("diffusions.jl")
export pagerank, pagerank_power!, personalized_pagerank, seeded_pagerank, stochastic_mult!


end # end module
