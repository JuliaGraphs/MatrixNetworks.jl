function mst_prim_test()
    A = load_matrix_network("airports")
    A = -A
    T = mst_prim_matrix(A)
    return true
end