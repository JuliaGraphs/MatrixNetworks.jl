using MatrixNetworks
using Test
using Random
using Statistics
using SparseArrays
using LinearAlgebra

# Todo
# 1. Add core MatrixNetworks tests



all_tests = [
             "matrixnetwork",
             "utility",
             "generators",
             "bfs",
             "biconnected",
             "bipartite_matching",
             "clustercoeffs",
             "corenums",
             "cosineknn",
             "csr_to_sparse",
             "dfs",
             "diffusions",
             "dijkstra",
             "dirclustercoeffs",
             "floydwarshall",
             "largest_component",
             "mst_prim",
             "network_formatting",
             "scomponents",
             "spectral",
             "sparse_to_csr",
             "triangles"
             ]

#for t in all_tests
for ti = 1:length(all_tests)
    t = all_tests[ti]
    test_name = join(["$(t)", "_test",".jl"])
    @show test_name
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
