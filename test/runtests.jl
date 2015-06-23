using MatrixNetworks
using Base.Test
using Compat

all_tests = ["dfs_test", "bfs_test", "bipartite_matching_test", "clustercoeffs_test", 
             "scomponents_test"]

for t in all_tests
    test_path = joinpath(Pkg.dir("MatrixNetworks"), "test", "$(t).jl")
    println("running $(test_path) ...")
    include(test_path)
end