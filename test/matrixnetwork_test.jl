function matrixnetwork_test()
    A = load_matrix_network("dfs_example")
    M = MatrixNetwork(A)
    B = sparse(M)
    @test A == B
    C = sparse_transpose(M)
    @test A' == C
    
    load_matrix_network_all("minnesota")
    load_matrix_network_all("U3A")
    
    #@show issym(sparse([0 1; 0 0]'))
    @test MatrixNetworks.is_undirected(MatrixNetwork(sparse([0 1; 1 0]))) == true
    @test MatrixNetworks.is_undirected(MatrixNetwork(sparse([0 1; 0 1]))) == false
    @test MatrixNetworks.is_undirected(MatrixNetwork(load_matrix_network("dfs_example"))) == false

    return true
    
    
end

