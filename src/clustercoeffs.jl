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

    M = sprand(A.n,A.n,0.0)
    M.colptr = A.rp
    M.rowval = A.ci
    M.nzval = A.vals
    
    if !(M' == M)
        error("Only undirected (symmetric) inputs are allowed")
    end
    
    (rp,ci,ai) = (A.rp,A.ci,A.vals)
    if length(find(ai.<0)) != 0
        error("only positive edge weights allowed")
    end
    return clustercoeffs_phase2(donorm,rp,ci,ai,usew)
end

function clustercoeffs_phase2{T}(donorm::Bool,rp::Vector{Int64},ci::Vector{Int64},
                                          ai::Vector{T},usew::Bool)
    n = length(rp) - 1
    cc = zeros(Float64,n)
    # ind = falses(n,1)
    ind = zeros(Bool,n)
    cache=zeros(Float64,n)
    ew = 1
    ew2 = 1
    for v = 1:n
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
        curcc = 0
        d = rp[v+1]-rp[v]
        # run two steps of bfs to try and find triangles. 
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            if v == w
                d = d-1
                continue
            end #discount self-loop
            for rpi2 = rp[w]:rp[w+1]-1
                x = ci[rpi2]
                if x == w
                    continue
                end
                if ind[x]
                    if usew
                        ew = ai[rpi]
                        ew2 = ai[rpi2]
                    end
                    curcc = curcc+ew^(1/3)*ew2^(1/3)*cache[x]
                end
            end
        end
        if donorm && d>1
            cc[v] = curcc/(d*(d-1))
        elseif d>1
            cc[v] = curcc
        end
        
        for rpi = rp[v]:rp[v+1]-1
            w = ci[rpi]
            ind[w] = 0
        end # reset indicator
    end
    return cc
end 

############################
### Additional functions ###
############################
# sparse matrices:
function clustercoeffs{T}(A::SparseMatrixCSC{T,Int64},weighted::Bool,normalized::Bool)
    donorm = true
    usew = true
    if !normalized
        donorm = false
    end
    if !weighted
        usew = false
    end
    
    At = A'
    if !(At == A)
        error("Only undirected (symmetric) inputs are allowed")
    end

    (rp,ci,ai) = (At.colptr,At.rowval,At.nzval);
    if length(find(ai.<0)) != 0
        error("only positive edge weights allowed")
    end
    return clustercoeffs_phase2(donorm,rp,ci,ai,usew)
end

function clustercoeffs{T}(A::SparseMatrixCSC{T,Int64})
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
function clustercoeffs{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64)
    return clustercoeffs(MatrixNetwork(n,rp,ci,vals))
end

function clustercoeffs{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64,weighted::Bool,normalized::Bool)
    return clustercoeffs(MatrixNetwork(n,rp,ci,vals),weighted,normalized)
end