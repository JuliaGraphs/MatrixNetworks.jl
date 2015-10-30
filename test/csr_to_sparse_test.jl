function csr_to_sparse_test()
    i = [1;2;3]
    j = [3;4;4]
    v = [8;9;10]
    (rp,ci,ai,m) = sparse_to_csr(i,j,v)
    (nzi,nzj,nzv) = csr_to_sparse(rp,ci,ai)
    A = sparse(nzi,nzj,nzv,length(rp)-1,maximum(ci))
    
    # more tests added here
    # clique to sparse test
    rp = 1:5:26
    rp = collect(rp)
    ci = vec(reshape(repmat(1:5,5,1)',25,1))
    ai = ones(Int64,25)
    (i,j,k) = csr_to_sparse(rp,ci,ai)
    if !isequal(full(M),ones(5,5))
        error("csr_to_sparse_test failed")
    end
    
    # 100 random trials
    



    return true
end