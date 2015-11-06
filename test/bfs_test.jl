function bfs_test()
    A = load_matrix_network("bfs_example")
    bfs(MatrixNetwork(A),1)
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 0]
    (d,dt,pred) = bfs(MatrixNetwork(sparse(A)),1)
    d_ref = [0,2,1,1]
    dt_ref = [0,3,1,2]
    pred_ref = [1,3,1,1]
    if d_ref != d_ref || dt_ref !=dt || pred_ref != pred
        error("bfs failed")
    end    
    return true
end