using MatrixNetworks
using Base.Test
using Compat
using Lint

all_tests = ["dfs", "bfs", "bipartite_matching", "clustercoeffs", 
             "scomponents","corenums","mst_prim"]

for t in all_tests
    test_path = joinpath(Pkg.dir("MatrixNetworks"), "test", join(["$(t)", "_test",".jl"])
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