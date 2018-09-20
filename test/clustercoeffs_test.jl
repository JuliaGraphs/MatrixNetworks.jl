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
        bad_mat = spdiagm(1=>[1;1]) #spdiagm(([1; 1],), 1, 3, 3)

        @test_throws ErrorException clustercoeffs(MatrixNetwork(bad_mat))
        @test_throws ErrorException clustercoeffs(bad_mat)
        # negative weights
        @test_throws ErrorException clustercoeffs(MatrixNetwork(-sparse(1.0I,4,4)))
        @test_throws ErrorException clustercoeffs(-sparse(1.0I,4,4))
    end
    cc_mn = clustercoeffs(MatrixNetwork(A), false, false)
    cc_csc = clustercoeffs(A, false, false)
    cc_csr = clustercoeffs(sparse_to_csr(A)..., false, false)
    cc_tri = clustercoeffs(findnz(A)[1], findnz(A)[2], false, false)
    @test cc_mn == cc_csc
    @test cc_csr == cc_csc
    @test cc_tri == cc_csc
end
