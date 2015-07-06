# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation
# TODO: add more examples
# TODO support MatrixNetworks input
# TODO support other similarity functions
# TODO allow option for diagonal similarity too
# TODO allow option for triplet output


"""
Example
-------

S = corenums(A)

cosineknn compute the k-nearest neighbors similarity metric between the
vertices of A or the upper half of a bipartite graph A
"""

function cosineknn{T}(A::SparseMatrixCSC{T,Int64},k::Int64)
    (rp,ci,ai)=sparse_to_csr(A)
    (rpt,cit,ait)=sparse_to_csr(A')
    (m,n) = size(A)
    # accumarray
    rn = zeros(Float64,maximum(cit))
    for i = 1:length(ait)
        rn[cit[i]] += ait[i]^2
    end
    rn = sqrt(rn)
    rn[rn.>eps()] = 1./rn[rn>eps(1)] # do division once
    
    si = zeros(Int64,m*K)
    si = zeros(Int64,m*K)
    sv = zeros(Int64,m*K)
    nused = 0
    
    dwork = zeros(Int64,m)
    iwork = zeros(Int64,m)
    iused = zeros(Int64,m)
    
    for i = 1:m
    
        # for each row of A, compute an inner product against all rows of A.
        # to do this efficiently, we know that the inner-product will
        # only be non-zero if two rows of A overlap in a column.  Thus, in
        # row i of A, we look at all non-zero columns, and then look at all
        # rows touched by those columns in the A' matrix.  We  compute
        # the similarity metric, then sort and only store the top k.
        
        curused = 0 # track how many non-zeros we compute during this operation
        for ri = rp[i]:rp[i+1]-1
            # find all columns j in row i
            j = ci[ri]
            aij = ai[ri]*rn[i]
            
            for rit = rpt[j]:rpt[j+1]-1
            
                # find all rows k in column j
                k = cit[rit]
                if k == i
                    continue
                end # skip over the standard entries
                akj = ait[rit]*rn[k]
                if iwork[k]>0
                    # we already have a non-zero between row i and row k
                    dwork[iwork[k] = dwork[iwork[k]] + aij*akj
                else
                    # we need to add a non-zero betwen row i and row k
                    curused = curused + 1
                    iwork[k] = curused
                    dwork[curused] = aij*akj
                    iused[curused] = k
                end
            end
        end
        
        # don't sort if fewer than K elements
        if curused < K
            sperm = 1:curused
        else
            sperm = sortperm(dwork[1:curused])
        end
    
        for k = 1:minimum(K,length(sperm))
            nused = nused + 1
            si[nused] = i
            sj[nused] = iused[sperm[k]]
            sv[nused] = dwork[sperm[k]]
        end
    
        # reset the iwork array, no need to reset dwork, as it's just 
        # overwritten
        for k = 1:curused
            iwork[iused[k]] = 0
        end
    end
    
    si = si[1:nused]
    sj = sj[1:nused]
    sv = sv[1:nused]
    
    S = sparse(si,sj,sv,n,n)
    
    return S
    
end