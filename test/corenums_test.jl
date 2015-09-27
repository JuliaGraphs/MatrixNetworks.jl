function corenums_test()
    A = load_matrix_network("cores_example")
    corenums(MatrixNetwork(A))
    return true
end