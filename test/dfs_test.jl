using MAT
function dfs_test()
    file_path = Pkg.dir("MatrixNetworks/data/dfs_example.mat")
    file = matopen(file_path)
    A = read(file,"A")
    close(file)
    return (dfs(MatrixNetwork(A),1))
end