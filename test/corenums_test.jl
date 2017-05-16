@testset "corenums" begin
    A = load_matrix_network("cores_example")
    A_ref = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 1]
    d_ref = [3;3;3;3]
    rt_ref = [1;2;3;4]
    @testset "MatrixNetwork" begin
        corenums(MatrixNetwork(A))
        (d,rt) = corenums(MatrixNetwork(sparse(A_ref)))
        @test d_ref == d
        @test rt_ref == rt
    end
    @testset "SparseMatrixCSC" begin
        corenums(sparse(A))
        (d,rt) = corenums(sparse(A_ref))
        @test d_ref == d
        @test rt_ref == rt
    end
    @testset "CSR" begin
        corenums(sparse_to_csr(sparse(A))...)
        (d,rt) = corenums(sparse_to_csr(sparse(A_ref))...)
        @test d_ref == d
        @test rt_ref == rt
    end
    @testset "triplet" begin
        corenums(findnz(sparse(A_ref))[1], findnz(sparse(A_ref))[2])
        (d,rt) = corenums(findnz(sparse(A_ref))[1], findnz(sparse(A_ref))[2])
        @test d_ref == d
        @test rt_ref == rt
    end
end
