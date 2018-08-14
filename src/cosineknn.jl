"""
COSINEKNN
---------
    compute the k-nearest neighbors similarity metric between the
    vertices of A or the upper half of a bipartite graph A.

Functions
---------
- S = cosineknn{T}(A::SparseMatrixCSC{T,Int64},K::Int64)
- S = cosineknn(A::MatrixNetwork,K::Int64)

Example
-------
~~~
A = load_matrix_network("bfs_example")
S = cosineknn(A,2)
~~~
"""
function cosineknn(A::SparseMatrixCSC{T,Int64},K::Int64) where T
    (rp,ci,ai) = sparse_to_csr(A)
    (rpt,cit,ait) = sparse_to_csr(copy(A'))
    (m,n) = size(A)
    return cosineknn_internal(rp,ci,ai,rpt,cit,ait,m,K)
end

##  MatrixNetwork support
function cosineknn(A::MatrixNetwork,K::Int64)
    # for the original matrix
    (rp,ci,ai) = (A.rp,A.ci,A.vals)
    # for the transposed matrix:
    M = sparse_transpose(A)
    At = MatrixNetwork(M)
    (rpt,cit,ait) = (At.rp,At.ci,At.vals)
    m = A.n
    return cosineknn_internal(rp,ci,ai,rpt,cit,ait,m,K)

end

function cosineknn_internal(rp::Vector{Int64},ci::Vector{Int64},ai::Vector{T},
                rpt::Vector{Int64},cit::Vector{Int64},ait::Vector{T},m::Int64,K::Int64) where T

    # accumarray
    rn = zeros(Float64,maximum(cit))
    for i = 1:length(ait)
        rn[cit[i]] += ait[i]^2
    end
    rn = sqrt.(rn)
    
     
    rn[rn.>eps()] = 1 ./rn[rn.>eps()] # do division once
    
    si = zeros(Int64,m*K)
    sj = zeros(Int64,m*K)
    sv = zeros(Float64,m*K)
    nused = 0
    
    dwork = zeros(Float64,m)
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
                    dwork[iwork[k]] += aij*akj
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
            sperm = sortperm(dwork[1:curused],rev=true)
        end
    
        for k = 1:min(K,length(sperm))
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
    
    ei = si[1:nused]
    ej = sj[1:nused]
    ev = sv[1:nused]
    
    S = sparse(ei,ej,ev,m,m)
    
    return S
end
