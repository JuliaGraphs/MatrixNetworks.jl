include("../src/manage_data.jl")
function bfs_test()
    file_path = Pkg.dir("MatrixNetworks/data/bfs_example.smat")
    A = readSMAT(file_path)
    return bfs(MatrixNetwork(A),1)
end