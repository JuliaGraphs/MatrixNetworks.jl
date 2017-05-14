function clustercoeffs_test()
    A = load_matrix_network("clique-10")
    cc = clustercoeffs(MatrixNetwork(A))

    v = ones(Int64,10)
    if v != cc
        error("clustercoeffs failed")
    end
   
    # SparseMatrixCSC
    cc = clustercoeffs(A)

    v = ones(Int64,10)
    if v != cc
        error("clustercoeffs failed")
    end
    
    # CSR 
    cc = clustercoeffs(sparse_to_csr(A)...)

    v = ones(Int64,10)
    if v != cc
        error("clustercoeffs failed")
    end
    
    # triplet
    cc = clustercoeffs(findnz(A)[1], findnz(A)[2])

    v = ones(Int64,10)
    if v != cc
        error("clustercoeffs failed")
    end

    # not undirected
    bad_mat = spdiagm(([1; 1],), 1, 3, 3)
    @test_throws ErrorException clustercoeffs(MatrixNetwork(bad_mat))
    @test_throws ErrorException clustercoeffs(bad_mat)
    # negative weights
    @test_throws ErrorException clustercoeffs(MatrixNetwork(-speye(4)))
    @test_throws ErrorException clustercoeffs(-speye(4))
    
    return true
end
