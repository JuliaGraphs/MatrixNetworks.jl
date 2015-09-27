function cosineknn_test()
    A = load_matrix_network("bfs_example")
    S = cosineknn(A,2)
    return true
end