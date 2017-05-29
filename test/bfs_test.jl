@testset "bfs" begin
    A = load_matrix_network("bfs_example")
    A_ref = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    
    @testset "MatrixNetwork" begin
        bfs(MatrixNetwork(A), 1)
        
        (d,dt,pred) = bfs(MatrixNetwork(sparse(A_ref)), 1)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    
        bfs(MatrixNetwork(A), 1, 0)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    end 
    @testset "SparseMatrixCSC" begin
        bfs(A,1)
    
        (d,dt,pred) = bfs(sparse(A_ref), 1)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    
        bfs(A,1,0)
        (d,dt,pred) = bfs(sparse(A_ref ),1,0)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    end 
    @testset "triplet" begin
        bfs(findnz(A)[1], findnz(A)[2], 1)
    
        (d,dt,pred) = bfs(findnz(sparse(A_ref))[1], findnz(sparse(A_ref))[2],1)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    
        bfs(findnz(A)[1], findnz(A)[2], 1, 0)
    
        (d,dt,pred) = bfs(findnz(sparse(A_ref))[1], findnz(sparse(A_ref))[2], 1, 0)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    end
    @testset "CSR" begin
        A_csr = sparse_to_csr(A)
        bfs(A_csr..., 1)
    
        (d,dt,pred) = bfs(sparse_to_csr(sparse(A_ref))..., 1)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    
        bfs(A_csr..., 1, 0)
    
        (d,dt,pred) = bfs(sparse_to_csr(sparse(A_ref))..., 1, 0)
        @test d == d_ref 
        @test dt_ref == dt
        @test pred_ref == pred
    end
end
