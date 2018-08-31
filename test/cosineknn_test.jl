@testset "cosineknn" begin
    A = load_matrix_network("bfs_example")
    S = cosineknn(A,2)
    
    # A = speye(Int64,4)
    A = sparse(1.0I,4,4)
    A[4,1] = 1
    CKNN = cosineknn(A,2)
    @test ((CKNN[4,1] *2)^2 - 2) <= 1e-7
    
    # A = speye(Int64,4)
    A = sparse(1.0I,4,4)
    A[4,1] = 1
    CKNN = cosineknn(MatrixNetwork(A),2)
    @test ((CKNN[4,1] *2)^2 - 2) <= 1e-7
end
