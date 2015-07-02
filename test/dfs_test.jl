include("../src/manage_data.jl")
function dfs_test()
    file_path = Pkg.dir("MatrixNetworks/data/dfs_example.smat")
    A = readSMAT(file_path)
    return (dfs(MatrixNetwork(A),1))
end