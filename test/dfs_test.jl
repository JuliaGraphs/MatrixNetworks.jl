function dfs_test()
    A = load_matrix_network("dfs_example")
    dfs(MatrixNetwork(A),1)
    
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = sparse(full(A))
    (d,dt,ft,pred) = dfs(B,1)
    if !(d == collect(0:n-1) && dt == collect(0:n-1) && ft == collect(2*n-1:-1:n) && 
                                    pred == collect(0:n-1)) 
        error("dfs test failed")
    end
    
    # sparse matrix
    A = load_matrix_network("dfs_example")
    dfs(A,1)
    
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = sparse(full(A))
    (d,dt,ft,pred) = dfs(B,1)
    if !(d == collect(0:n-1) && dt == collect(0:n-1) && ft == collect(2*n-1:-1:n) && 
                                    pred == collect(0:n-1)) 
        error("dfs test failed")
    end
    
    # sparse csr matrix
    A = load_matrix_network("dfs_example")
    dfs(sparse_to_csr(A)..., 1)
    
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = sparse(full(A))
    (d,dt,ft,pred) = dfs(B,1)
    if !(d == collect(0:n-1) && dt == collect(0:n-1) && ft == collect(2*n-1:-1:n) && 
                                    pred == collect(0:n-1)) 
        error("dfs test failed")
    end
    
    # triplet 
    A = load_matrix_network("dfs_example")
    dfs(findnz(A)[1], findnz(A)[2], 1)
    
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = sparse(full(A))
    (d,dt,ft,pred) = dfs(B,1)
    if !(d == collect(0:n-1) && dt == collect(0:n-1) && ft == collect(2*n-1:-1:n) && 
                                    pred == collect(0:n-1)) 
        error("dfs test failed")
    end
    return true
end
