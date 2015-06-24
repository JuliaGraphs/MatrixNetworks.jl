module MatrixNetworks

using Compat

if VERSION < v"0.4-"
    using Docile
end

""" 
Module ``MatrixNetworks``: Documentation on the module 

- Option 1: start with a sparse matrix A:
- example: ``M = MatrixNetwork(A)``
- Option 2: start with row and column indexes for nonzeros in the matrix
- example: ``M = MatrixNetwork(ei,ej)``
"""
MatrixNetworks

include("MatrixNetwork.jl")

function MatrixNetwork(A::SparseMatrixCSC{Float64,Int64})
    At = A'
    return MatrixNetwork(size(At,2),At.colptr,At.rowval,At.nzval)
end

function MatrixNetwork(ei::Vector{Int64},ej::Vector{Int64}) 
    At = sparse(ej,ei,true);
    return MatrixNetwork(size(At,2),At.colptr,At.rowval,At.nzval)
end

include("scomponents.jl")
include("bipartite_matching.jl")
include("bfs.jl")
include("dfs.jl")
include("clustercoeffs.jl")
include("corenums.jl")
include("readSMAT.jl")

# export everything to make them accessible as functions
export MatrixNetwork, bipartite_matching, bfs, dfs, clustercoeffs, 
corenums, scomponents, strong_components_map, readSMAT, enrich

# examples:
# include("../test/bfs_test.jl")
# include("../test/bipartite_matching_test.jl")
# include("../test/clustercoeffs_test.jl")
# include("../test/corenums_test.jl")
# include("../test/dfs_test.jl")
# include("../test/scomponents_test.jl")

# export test cases?
# export bfs_test

end # end module
