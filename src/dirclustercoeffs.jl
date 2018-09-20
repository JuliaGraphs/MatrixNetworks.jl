"""
DIRCLUSTERCOEFFS
----------------
    compute clustering coefficients for a directed graph.
    cc = dirclustercoeffs(A) returns the directed clustering coefficients
    (which generalize the clustering coefficients of an undirected graph, 
    and so calling this function on an undirected graph will produce the same
    answer as clustercoeffs, but less efficiently.)
    
    This function implements the algorithm from Fagiolo, Phys Rev. E. 
    76026107 (doi:10:1103/PhysRevE.76.026107).  
    
    (cc,cccyc,ccmid,ccin,ccout,nf) = dirclusteringcoeffs(A) returns different 
    components of the clustering coefficients corresponding to cycles,
    middles, in triangles, and out triangles.  See the manuscript for a 
    description of the various types of triangles counted in the above metrics.

Functions
---------
- (cc,cccyc,ccmid,ccin,ccout,nf) = dirclustercoeffs{T}(A::SparseMatrixCSC{T,Int64},weighted::Bool)
- (cc,cccyc,ccmid,ccin,ccout,nf) = dirclustercoeffs{T}(A::SparseMatrixCSC{T,Int64})
- (cc,cccyc,ccmid,ccin,ccout,nf) = dirclustercoeffs{T}(A::SparseMatrixCSC{T,Int64},weighted::Bool,normalized::Bool)

Example
-------
~~~
(A,xy,labels) = load_matrix_network_metadata("celegans")
(cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A, true, true)
(maxval, maxind) = findmax(cc)
labels[maxind]
~~~
"""
function dirclustercoeffs(A::SparseMatrixCSC{T,Int64},weighted::Bool) where T
    return dirclustercoeffs(A,weighted,true)
end

function dirclustercoeffs(A::SparseMatrixCSC{T,Int64}) where T
    return dirclustercoeffs(A,true,true)
end

function dirclustercoeffs(A::SparseMatrixCSC{T,Int64},weighted::Bool,normalized::Bool) where T
    donorm = true
    usew = true
    
    if !normalized
        donorm = false
    end
    if !weighted
        usew = false
    end
    
    if usew
        (rp,ci,ai) = sparse_to_csr(A)
        (col_ptr,ri,ati) = sparse_to_csr(copy(A'))
    else
        (rp,ci) = sparse_to_csr(A)
        (col_ptr,ri) = sparse_to_csr(copy(A'))
    end
    
    if any(ai.<0)
        error("only positive edge weights allowed")
    end

    n = length(rp)-1
    # initialize all the variables
    cc = zeros(Float64,n)
    ind = zeros(Bool,n)
    cache = zeros(Float64,n)
    degs = zeros(Int64,n)
    cccyc = zeros(Float64,n)
    ccmid = zeros(Float64,n)
    ccin = zeros(Float64,n)
    ccout = zeros(Float64,n)
    nf = zeros(Int64,n)
    
    # precompute degrees
    for v = 1:n
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v == w
                continue
            else
                degs[w] = degs[w] + 1
                degs[v] = degs[v] + 1
            end
        end
    end
    
    ew = one(T)
    ew2 = one(T)
    for v = 1:n
        # setup counts for the different cycle types
        bilatedges = 0.0
        curcccyc = 0.0
        curccmid = 0.0
        curccin = 0.0
        curccout = 0.0
        # 1.  
        # find triangles with out links as last step, so precompute the inlinks
        # back to node v
        for cpi = col_ptr[v]:col_ptr[v+1]-1
            w = ri[cpi]
            if usew
                ew = ati[cpi]
            end
            if v != w
                ind[w] = 1
                cache[w] = ew^(1/3)
            end
        end
        # count cycles (cycles are out->out->out)
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v == w
                continue
            end # discount self-loop
            if usew
                ew = ai[rpi]^(1/3)
            end
            for rpi2 = rp[w]:rp[w+1]-1
                x = ci[rpi2]
                if x == w
                    continue
                end
                if x == v
                    bilatedges = bilatedges + 1
                    continue
                end
                if ind[x] == 1
                    if usew
                        ew2 = ai[rpi2]
                    end
                    curcccyc = curcccyc + ew*ew2^(1/3) * cache[x]
                end
            end
        end
        # count middle-man circuits (out->in->out)
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v == w
                continue
            end # discount self-loop
            if usew
                ew = ai[rpi]^(1/3)
            end
            for cpi = col_ptr[w]:col_ptr[w+1]-1
                x=ri[cpi]
                if x==w
                    continue
                end
                if ind[x] == 1
                    if usew
                        ew2 = ati[cpi]
                    end
                    curccmid=curccmid+ew*ew2^(1/3)*cache[x]
                end
            end
        end
        # count in-link circuits (in->out->out)
        for cpi = col_ptr[v]:col_ptr[v+1]-1
            w = ri[cpi]
            if v == w
                continue
            end # discount self-loop
            if usew
                ew = ati[cpi]^(1/3)
            end
            for rpi2 = rp[w]:rp[w+1]-1
                x = ci[rpi2]
                if x == w
                    continue
                end
                if ind[x] == 1
                    if usew
                        ew2 = ai[rpi2]
                    end
                    curccin = curccin + ew*ew2^(1/3) * cache[x]
                end
            end
        end
        # reset and reinit the cache for outlinks
        for cpi = col_ptr[v]:col_ptr[v+1]-1
            w = ri[cpi]
            ind[w]=0
        end # reset indicator
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if usew
                ew = ai[rpi]
            end
            if v != w
                ind[w] = 1
                cache[w] = ew^(1/3)
            end
        end
        # count out-link circuits (out->out->in)
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v == w
                continue
            end # discount self-loop
            if usew
                ew = ai[rpi]^(1/3)
            end
            for rpi2 = rp[w]:rp[w+1]-1
                x = ci[rpi2]
                if x == w
                    continue
                end
                if ind[x] == 1
                    if usew
                        ew2 = ai[rpi2]
                    end
                    curccout = curccout+ew*ew2^(1/3) * cache[x]
                end
            end
        end
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            ind[w] = 0
        end # reset indicator
        # store the values
        curnf = degs[v]*(degs[v]-1) - 2*bilatedges
        curcc = curcccyc + curccmid + curccin + curccout
        if curnf>0 && donorm
            curcc = curcc/curnf
        end
        cc[v] = curcc
        cccyc[v] = curcccyc
        ccmid[v] = curccmid
        ccin[v] = curccin
        ccout[v] = curccout
        nf[v] = curnf

    end
    return (cc,cccyc,ccmid,ccin,ccout,nf)
end
