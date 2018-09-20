@testset "dfs" begin
    A = load_matrix_network("dfs_example")
    dfs(MatrixNetwork(A),1)
    
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A_ref = Tridiagonal(O,Z,O)
    B = sparse(Matrix(A_ref))
    @testset "SparseMatrixCSC" begin
        dfs(A,1)
        
        (d,dt,ft,pred) = dfs(B,1)
        @test d == collect(0:n-1)
        @test dt == collect(0:n-1)
        @test ft == collect(2*n-1:-1:n)
        @test pred == collect(0:n-1)
    end 
    @testset "CSR" begin
        dfs(sparse_to_csr(A)..., 1)
    
        (d,dt,ft,pred) = dfs(sparse_to_csr(B)..., 1)
        @test d == collect(0:n-1)
        @test dt == collect(0:n-1)
        @test ft == collect(2*n-1:-1:n)
        @test pred == collect(0:n-1)
    end 
    @testset "triplet" begin
        dfs(findnz(A)[1], findnz(A)[2], 1)
        
        (d,dt,ft,pred) = dfs(findnz(B)[1], findnz(B)[2], 1)
        @test d == collect(0:n-1)
        @test dt == collect(0:n-1)
        @test ft == collect(2*n-1:-1:n)
        @test pred == collect(0:n-1)
    end
end
