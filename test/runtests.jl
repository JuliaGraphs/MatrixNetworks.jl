using MatrixNetworks
using Base.Test
#using Lint

# Todo
# 1. Add core MatrixNetworks tests



all_tests = ["matrixnetwork",
             "generators",
             "bfs",
             "biconnected",
             "ear_decomposition",
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
             "scomponents",
             "spectral",
             "sparse_to_csr"]

for t in all_tests
    test_name = join(["$(t)", "_test",".jl"])
    test_path = joinpath(dirname(@__FILE__), test_name)
    println("running $(test_path) ...")
    test_function = include(test_path)
    test_function() 
    println("running $(test_path) ... PASSED")
end



#println("testing package with Lint...")
#msgs = lintpkg( "MatrixNetworks", returnMsgs = true )
#if isempty(msgs)
#    info("Lint package passed")
#else
#    warn("Lint package didn't pass")
#end
