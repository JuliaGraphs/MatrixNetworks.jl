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

function MatrixNetwork{T}(A::SparseMatrixCSC{T,Int64})
    At = A'
    return MatrixNetwork(size(At,2),At.colptr,At.rowval,At.nzval)
end

function MatrixNetwork(ei::Vector{Int64},ej::Vector{Int64}) 
    At = sparse(ej,ei,true);
    return MatrixNetwork(size(At,2),At.colptr,At.rowval,At.nzval)
end

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
include("largest_component.jl")
include("cosineknn.jl")
include("dirclustercoeffs.jl")
include("dijkstra.jl")
include("mst_prim.jl")

# export everything to make them accessible as functions
export MatrixNetwork, bipartite_matching, edge_list, create_sparse,
bipartite_matching_setup, bipartite_matching_indicator, bfs, 
dfs, clustercoeffs, corenums, scomponents, strong_components_map, 
readSMAT, enrich, load_matrix_network, matrix_network_datasets, 
csr_to_sparse, load_matrix_network_metadata, floydwarshall, largest_component,
sparse_to_csr, cosineknn, dirclustercoeffs, dijkstra, mst_prim, mst_prim_matrix,
csr_to_sparse_matrix


end # end module
