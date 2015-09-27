function largest_component_test()
    A = load_matrix_network("dfs_example")
    (Acc,p) = largest_component(A)
    return true
end