function floydwarshall_test()
    A = load_matrix_network("all_shortest_paths_example")
    (D,P) = floydwarshall(A)
    return true
end