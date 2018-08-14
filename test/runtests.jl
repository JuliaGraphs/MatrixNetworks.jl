using MatrixNetworks
using Test
using Random
using Statistics
using SparseArrays
using Arpack
#using Lint

# Todo
# 1. Add core MatrixNetworks tests



all_tests = [
             #"matrixnetwork",
             #"generators",
             #"bfs",
             #"biconnected",
             #"bipartite_matching",
             #"clustercoeffs",
             #"corenums",
             #"cosineknn",
             #"csr_to_sparse",
             #"dfs",
             "diffusions", # not tested yet
             #"dijkstra",
             #"dirclustercoeffs",
             #"floydwarshall",
             #"largest_component",
             #"mst_prim",
             #"scomponents",
             #"spectral",
             #"sparse_to_csr"
             ]

for t in all_tests
    test_name = join(["$(t)", "_test",".jl"])
    test_path = joinpath(dirname(@__FILE__), test_name)
    include(test_path)
end



#println("testing package with Lint...")
#msgs = lintpkg( "MatrixNetworks", returnMsgs = true )
#if isempty(msgs)
#    info("Lint package passed")
#else
#    warn("Lint package didn't pass")
#end
