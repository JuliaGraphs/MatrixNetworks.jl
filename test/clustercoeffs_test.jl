function clustercoeffs_test()
    A = load_matrix_network("clique-10")
    clustercoeffs(MatrixNetwork(A))
    return true
end