@testset "clustercoeffs" begin
    A = load_matrix_network("clique-10")
    v = ones(Int64,10)
    @testset "MatrixNetwork" begin
        cc = clustercoeffs(MatrixNetwork(A))
        @test v == cc
    end
    @testset "SparseMatrixCSC" begin
        cc = clustercoeffs(A)
        @test v == cc
    end
    @testset "CSR" begin 
        cc = clustercoeffs(sparse_to_csr(A)...)
        @test v == cc
    end
    @testset "triplet" begin
        cc = clustercoeffs(findnz(A)[1], findnz(A)[2])
        @test v == cc
    end
    @testset "error throwing" begin
        # not undirected
        bad_mat = spdiagm(([1; 1],), 1, 3, 3)
        @test_throws ErrorException clustercoeffs(MatrixNetwork(bad_mat))
        @test_throws ErrorException clustercoeffs(bad_mat)
        # negative weights
        @test_throws ErrorException clustercoeffs(MatrixNetwork(-speye(4)))
        @test_throws ErrorException clustercoeffs(-speye(4))
    end
end
