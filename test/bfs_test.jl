function bfs_test()
    A = load_matrix_network("bfs_example")
    bfs(MatrixNetwork(A), 1)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(MatrixNetwork(sparse(A)), 1)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    
    A = load_matrix_network("bfs_example")
    bfs(MatrixNetwork(A), 1, 0)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(MatrixNetwork(sparse(A)), 1, 0)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    # test for SparseMatrixCSC
    A = load_matrix_network("bfs_example")
    bfs(A,1)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(sparse(A), 1)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    
    A = load_matrix_network("bfs_example")
    bfs(A,1,0)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(sparse(A),1,0)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    
    # test for triplet
    A = load_matrix_network("bfs_example")
    bfs(findnz(A)[1], findnz(A)[2], 1)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(findnz(sparse(A))[1], findnz(sparse(A))[2],1)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    
    A = load_matrix_network("bfs_example")
    bfs(findnz(A)[1], findnz(A)[2], 1, 0)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(findnz(sparse(A))[1], findnz(sparse(A))[2], 1, 0)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end

    # test for CSR
    A = load_matrix_network("bfs_example")
    A_csr = sparse_to_csr(A)
    bfs(A_csr..., 1)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(sparse_to_csr(sparse(A))..., 1)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    
    A = load_matrix_network("bfs_example")
    A_csr = sparse_to_csr(A)
    bfs(A_csr..., 1, 0)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(sparse_to_csr(sparse(A))..., 1, 0)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end
    return true
end
