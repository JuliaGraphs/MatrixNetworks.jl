using MatrixNetworks
using Base.Test
using Compat
using Lint

all_tests = ["bfs",
             "bipartite_matching",
             "clustercoeffs",
             "corenums",
             "cosineknn",
             "csr_to_sparse",
             "dfs",
             "dijkstra",
             "dirclustercoeffs",
             "floydwarshall",
             "largest_component",
             "mst_prim",
             "scomponents",
             "sparse_to_csr"]

for t in all_tests
    test_name = join(["$(t)", "_test",".jl"])
    test_path = joinpath(Pkg.dir("MatrixNetworks"), "test", test_name)
    println("running $(test_path) ...")
    test_function = include(test_path)
    @test test_function() == true
end

println("testing package with Lint...")
msgs = lintpkg( "MatrixNetworks", returnMsgs = true )
if isempty(msgs)
    info("Lint package passed")
else
    warn("Lint package didn't pass")
end