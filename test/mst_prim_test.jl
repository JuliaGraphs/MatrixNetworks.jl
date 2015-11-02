function mst_prim_test()
    A = load_matrix_network("airports")
    A = -A
    T = mst_prim_matrix(A)
    
    A = load_matrix_network("clr-24-1")
    A[2,3] = 9
    A[3,2] = 9
    T = mst_prim_matrix(A)
    Ttrueijv = [
     2     8     1     4     6     9     3     5     4     3     7     6     8     1     7     3
     1     1     2     3     3     3     4     4     5     6     6     7     7     8     8     9
     4     8     4     7     4     2     7     9     9     4     2     2     1     8     1     2 ]
    Ttrue = sparse(vec(Ttrueijv[1,:]),  vec(Ttrueijv[2,:]), vec(Ttrueijv[3,:]), 9,9)
     
    if nnz(T - Ttrue) != 0
        error("mst_prim failed")
    end
    
    A = load_matrix_network("clr-24-1")
    T1 = mst_prim_matrix(A)
    T2 = mst_prim_matrix(A,false,5)
    T1T2diff = sparse([2,1],[3,8],[8,-8],9,9)
    
    if nnz(triu(T1-T2) - T1T2diff) != 0
        error("mst_prim failed")
    end

    return true
end