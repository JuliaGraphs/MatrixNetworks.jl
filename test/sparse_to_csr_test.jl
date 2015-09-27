function sparse_to_csr_test()
    i = [1;2;3]
    j = [3;4;4]
    v = [8;9;10]
    (rp,ci,ai,m) = sparse_to_csr(i,j,v,4)
    return true
end