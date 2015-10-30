function sparse_to_csr_test()
    i = [1;2;3]
    j = [3;4;4]
    v = [8;9;10]
    (rp,ci,ai,m) = sparse_to_csr(i,j,v)
    
    # more tests added here
    for t = 1:100
        A = sprand(100,80,0.01)
        (rp,ci,ai,m)=sparse_to_csr(A)
        i = zeros(Int64,length(ai))
        j = zeros(Int64,length(ai))
        a = zeros(Float64,length(ai))
        n = length(rp)-1
        nz = 0
        for cr = 1:n
            for ri = rp[cr]:rp[cr+1]-1
                nz=nz+1
                i[nz]=cr
                j[nz]=ci[ri]
                a[nz]=ai[ri]
            end
        end
        A2 = sparse(i,j,a,n,m)
        if ~isequal(A,A2)
            error("sparse_to_csr_test failed")
        end
    end
    return true
end
