function scomponents_test()
    A = load_matrix_network("cores_example")
    scomponents(MatrixNetwork(A))
    return true
end