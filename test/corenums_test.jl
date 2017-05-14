function corenums_test()
    A = load_matrix_network("cores_example")
    corenums(MatrixNetwork(A))
    
    A = [1 0 1 1; 0 1 1 1; 1 1 0 1; 1 1 1 1]
    A = sparse(A)
    (d,rt) = corenums(A)
    d_ref = [3;3;3;3]
    rt_ref = [1;2;3;4]
    if d_ref != d || rt_ref != rt
        error("corenums failed")
    end
    
    (d,rt) = corenums(sparse_to_csr(A)...)
    d_ref = [3;3;3;3]
    rt_ref = [1;2;3;4]
    if d_ref != d || rt_ref != rt
        error("corenums failed")
    end
    
    (d,rt) = corenums(findnz(A)[1], findnz(A)[2])
    d_ref = [3;3;3;3]
    rt_ref = [1;2;3;4]
    if d_ref != d || rt_ref != rt
        error("corenums failed")
    end
    return true
end
