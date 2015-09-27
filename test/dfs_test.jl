function dfs_test()
    A = load_matrix_network("dfs_example")
    dfs(MatrixNetwork(A),1)
    return true
end