function clustercoeffs_test()
    A = load_matrix_network("clique-10")
    cc = clustercoeffs(MatrixNetwork(A))
    
    v = ones(Int64,10)
    if v != cc
        error("clustercoeffs failed")
    end
    
    return true
end