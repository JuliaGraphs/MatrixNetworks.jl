using MAT
function clustercoeffs_test()
    file_path = Pkg.dir("MatrixNetworks/data/clique-10.mat")
    file = matopen(file_path)
    A = read(file,"A")
    close(file)
    cc = clustercoeffs(MatrixNetwork(A))
    return cc;
end