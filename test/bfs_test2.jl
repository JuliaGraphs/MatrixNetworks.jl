using MAT
function bfs_test2()
    file_path = Pkg.dir("MatrixNetworks/data/bfs_example.mat")
    file = matopen(file_path)
    A = read(file,"A")
    close(file)
    return bfs(MatrixNetwork(A),1)
end