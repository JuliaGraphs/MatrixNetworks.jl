function ear_test()

    A = load_matrix_network("ear_example")
    B = MatrixNetwork(A)
    count = ear_decomposition(B).non_trivial_ears
    @test count == 2

    A = load_matrix_network("clique-10")
    B = MatrixNetwork(A)
    count = ear_decomposition(B).non_trivial_ears
    @test count == 1

    A = empty_graph()
    count = ear_decomposition(A).non_trivial_ears
    @test count == 0

    return true
end
