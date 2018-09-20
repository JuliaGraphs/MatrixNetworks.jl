"""
CLUSTERCOEFFS
-------------
    compute undirected clustering coefficients for a graph. clustercoeffs(A) computes a 
    normalized, weighted clustering coefficients from a graph represented by a symmetric 
    adjacency matrix A. clustercoeffs(A,weighted,normalized), with weighted and normalized 
    boolean values indicating whether the computation has to be weighted and/or normalized.

Functions
---------
- cc = clustercoeffs(A::MatrixNetwork,weighted::Bool,normalized::Bool)
- cc = clustercoeffs{T}(A::SparseMatrixCSC{T,Int64},weighted::Bool,normalized::Bool)\n
If weighted and normalized are not specified, they are understood as true

Example
-------
~~~
A = load_matrix_network("clique-10")    
cc = clustercoeffs(MatrixNetwork(A))    
~~~
"""
function clustercoeffs(A::MatrixNetwork)
    return clustercoeffs(A, true, true);
end

function clustercoeffs(A::MatrixNetwork,weighted::Bool,normalized::Bool)
    donorm = true
    usew = true
    if !normalized
        donorm = false
    end
    if !weighted
        usew = false
    end

    if is_undirected(A) == false
        error("Only undirected (symmetric) inputs are allowed")
    end
    
    (rp,ci,ai) = (A.rp,A.ci,A.vals)
    if typeof(findfirst(ai.<0)) != Nothing
        error("only positive edge weights allowed")
    end
    return clustercoeffs_phase2(donorm,rp,ci,ai,usew)
end

function clustercoeffs_phase2(donorm::Bool,rp::Vector{Int64},ci::Vector{Int64},
                                       ai::Vector{T}, usew::Bool) where T
    n = length(rp) - 1
    cc = Vector{Float64}(undef,n)
    ind = zeros(Bool,n)
    cache = zeros(Float64,usew ? n : 0)

    @inbounds for v = 1:n
        
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v != w
                ind[w] = 1
                if usew # as of 2016-10-04, this makes a/insignificant slowdown
                    cache[w] = ai[rpi]^(1/3.0)
                end
            end
        end
        curcc = 0.0
        d = rp[v+1]-rp[v]
        # run two steps of bfs to try and find triangles. 
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v == w
                d = d-1
                continue
            end #discount self-loop
            
            istart=rp[w]
            iend = rp[w+1]-1
            if usew # as of 2016-10-04, this arrangement with outer if was better
                for rpi2 = istart:iend
                    x = ci[rpi2]
                    if ind[x]
                        curcc += ai[rpi]^(1/3)*ai[rpi2]^(1/3)*cache[x]*(x != w)
                    end
                end
            else
                for rpi2 = istart:iend
                    x = ci[rpi2]
                    if ind[x] # 
                        curcc += 1.0*(x != w)
                    end
                end
            end                 
            #=
            for rpi2 = rp[w]:rp[w+1]-1
                x = ci[rpi2]
                #if x == w
                    #continue
                #end
                if ind[x]
                    #if usew
                        #curcc += ai[rpi]^(1/3)*ai[rpi2]^(1/3)*cache[x]
                    #else
                    curcc += 1.0*(x != w)
                    #end
                end
            end
            =#
        end
        if donorm && d>1
            cc[v] = curcc/(d*(d-1))
        elseif d>1
            cc[v] = curcc
        end
        
        for rpi = rp[v]:rp[v+1]-1
            ind[ci[rpi]] = 0
        end # reset indicator
    end
    return cc
end 

############################
### Additional functions ###
############################
# sparse matrices:
function clustercoeffs(A::SparseMatrixCSC{T,Int64},weighted::Bool,normalized::Bool) where T
    donorm = true
    usew = true
    if !normalized
        donorm = false
    end
    if !weighted
        usew = false
    end
    
    if !is_undirected(A)
        error("Only undirected (symmetric) inputs are allowed")
    end

    (rp,ci,ai) = (A.colptr,A.rowval,A.nzval)
    if typeof(findfirst(ai.<0)) != Nothing
        error("only positive edge weights allowed")
    end
    return clustercoeffs_phase2(donorm,rp,ci,ai,usew)
end

function clustercoeffs(A::SparseMatrixCSC{T,Int64}) where T
    return clustercoeffs(A, true, true);
end

## Triplet Format:
function clustercoeffs(ei::Vector{Int64},ej::Vector{Int64})
    return clustercoeffs(MatrixNetwork(ei,ej))
end

function clustercoeffs(ei::Vector{Int64},ej::Vector{Int64},weighted::Bool,normalized::Bool)
    return clustercoeffs(MatrixNetwork(ei,ej),weighted,normalized)
end

## CSR sparse matrices - basically just like type MatrixNetwork
function clustercoeffs(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64) where T
    return clustercoeffs(MatrixNetwork(n,rp,ci,vals))
end

function clustercoeffs(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64,weighted::Bool,normalized::Bool) where T
    return clustercoeffs(MatrixNetwork(n,rp,ci,vals),weighted,normalized)
end