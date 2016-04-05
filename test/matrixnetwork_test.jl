function matrixnetwork_test()
    A = load_matrix_network("dfs_example")
    M = MatrixNetwork(A)
    B = sparse(M)
    @test A == B
    C = sparse_transpose(M)
    @test A' == C
    
    load_matrix_network_metadata_all("minnesota")
    load_matrix_network_metadata_all("U3A")

    return true
end

