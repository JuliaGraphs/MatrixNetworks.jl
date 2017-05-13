function cosineknn_test()
    A = load_matrix_network("bfs_example")
    S = cosineknn(A,2)
    
    A = speye(Int64,4)
    A[4,1] = 1
    CKNN = cosineknn(A,2)
    if ((CKNN[4,1] *2)^2 - 2) > 1e-7
        error("cosine knn failed")
    end
    
    A = speye(Int64,4)
    A[4,1] = 1
    CKNN = cosineknn(MatrixNetwork(A),2)
    if ((CKNN[4,1] *2)^2 - 2) > 1e-7
        error("cosine knn failed")
    end
    return true
end
