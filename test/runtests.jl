using MatrixNetworks
using Base.Test
using Compat
using Lint

all_tests = ["dfs_test", "bfs_test", "bipartite_matching_test", "clustercoeffs_test", 
             "scomponents_test","corenums_test"]

for t in all_tests
    test_path = joinpath(Pkg.dir("MatrixNetworks"), "test", "$(t).jl")
    println("running $(test_path) ...")
    include(test_path)
end

println("testing package with Lint...")
msgs = lintpkg( "MatrixNetworks", returnMsgs = true )
if isempty(msgs)
    println("Lint package passed")
else
    warn("Lint package didn't pass")
end