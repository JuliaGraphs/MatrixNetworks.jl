function csr_to_sparse_test()
    i = [1;2;3]
    j = [3;4;4]
    v = [8;9;10]
    (rp,ci,ai,m) = sparse_to_csr(i,j,v)
    (nzi,nzj,nzv) = csr_to_sparse(rp,ci,ai)
    A = sparse(nzi,nzj,nzv,length(rp)-1,maximum(ci))
end