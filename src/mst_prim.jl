"""
MST_PRIM
--------
    compute a minimum spanning tree with Prim's algorithm.
    T = mst_prim_matrix(A) computes a minimum spanning tree T using Prim's algorithm
    for the spanning tree of a graph with non-negative edge weights.
    
    T = mst_prim_matrix(A,false) produces an MST for just the component at A containing
    vertex 1.  T = mst_prim_matrix(A,0,u) produces the MST for the component
    containing vertex u.
    
    (ti,tj,tv,nverts) = mst_prim(A) returns the edges from the matrix and does not
    convert to a sparse matrix structure.  This saves a bit of work and is
    required when there are 0 edge weights.

Functions
---------
- (ti,tj,tv,nverts) = mst_prim(A::MatrixNetwork,full::Bool,u::Int64)
- (ti,tj,tv,nverts) = mst_prim(A::MatrixNetwork)
- (ti,tj,tv,nverts) = mst_prim{T}(A::SparseMatrixCSC{T,Int64},full::Bool,u::Int64)
- (ti,tj,tv,nverts) = mst_prim{T}(A::SparseMatrixCSC{T,Int64})
- T = mst_prim_matrix(A::MatrixNetwork,full::Bool,u::Int64) #output modifier
- T = mst_prim_matrix(A::MatrixNetwork) #output modifier
- T = mst_prim_matrix{T}(A::SparseMatrixCSC{T,Int64},full::Bool,u::Int64) #output modifier
- T = mst_prim_matrix{T}(A::SparseMatrixCSC{T,Int64}) #output modifier

Example
-------
~~~
A = load_matrix_network("airports")
A = -A #convert to travel time
A = max.(A,A')
A = sparse(A)
(ti,tj,tv,nverts) = mst_prim(A)
~~~
"""
function mst_prim(A::MatrixNetwork,full::Bool,u::Int64)
    
    (rp,ci,ai) = (A.rp,A.ci,A.vals)
    if any(ai.<0)
        error("prim''s algorithm cannot handle negative edge weights.")
    end
    nverts = length(rp) - 1
    d = Inf*ones(Float64,nverts)
    T = zeros(Int64,nverts)
    L = zeros(Int64,nverts)
    pred = zeros(Int64,length(rp)-1)
    
    # enter the main dijkstra loop
    for iter = 1:nverts
        if iter == 1
            root = u
        else
            root = mod(u+iter-1,nverts) + 1
            if L[root] > 0
                continue
            end
        end
        n = 1
        T[n] = root
        L[root] = n # oops, n is now the size of the heap
        d[root] = 0
        while n > 0
            v = T[1]
            L[v] = -1
            ntop = T[n]
            T[1] = ntop
            n = n-1
            if n > 0
                L[ntop] = 1
            end         # pop the head off the heap
            k = 1
            kt = ntop     # move element T[1] down the heap
            while true
                i = 2*k
                if i > n
                    break
                end       # end of heap
                if i == n
                    it = T[i]     # only one child, so skip
                else              # pick the smallest child
                    lc = T[i]
                    rc = T[i+1]
                    it = lc
                    if d[rc] < d[lc]
                        i = i+1
                        it = rc
                    end # right child is smaller
                end
                if d[kt] < d[it]
                    break     # at correct place, so end
                else
                    T[k] = it
                    L[it] = k
                    T[i] = kt
                    L[kt] = i
                    k=i # swap
                end
            end       # end heap down
            
            # for each vertex adjacent to v, relax it
            for ei = rp[v]:rp[v+1]-1      # ei is the edge index
                w = ci[ei]
                ew = ai[ei]          # w is the target, ew is the edge weight
                if L[w] < 0
                    continue
                end      # make sure we don't visit w twice
                # relax edge (v,w,ew)
                if d[w] > ew
                    d[w] = ew
                    pred[w] = v
                    # check if w is in the heap
                    k = L[w]
                    onlyup = false
                    if k == 0
                        # element not in heap, only move the element up the heap
                        n = n + 1
                        T[n] = w
                        L[w] = n
                        k = n
                        kt = w
                        onlyup = true
                    else
                        kt = T[k]
                    end
                    # update the heap, move the element down in the heap
                    while !onlyup
                        i = 2*k
                        if i > n
                            break
                        end          # end of heap
                        if i == n
                            it = T[i]    # only one child, so skip
                        else             # pick the smallest child
                            lc = T[i]
                            rc = T[i+1]
                            it = lc
                            if d[rc] < d[lc]
                                i = i+1
                                it = rc
                                end # right child is smaller
                            end
                            if d[kt] < d[it]
                                break    # at correct place, so end
                            else
                                T[k] = it
                                L[it] = k
                                T[i] = kt
                                L[kt] = i
                                k = i # swap
                            end
                        end
                        # move the element up the heap
                        j = k
                        tj = T[j]
                        while j > 1                      # j==1 => element at top of heap
                            j2 = convert(Int64,floor(j/2))
                            tj2 = T[j2]                  # parent element
                            if d[tj2] < d[tj]
                                break                    # parent is smaller, so done
                            else                         # parent is larger, so swap
                            T[j2] = tj
                            L[tj] = j2
                            T[j] = tj2
                            L[tj2] = j
                            j = j2
                        end
                    end
                end
            end
        end
        if !full
            break
        end
    end
    nmstedges = 0
    for i = 1:nverts
        if pred[i] > 0
            nmstedges = nmstedges + 1
        end
    end
    ti = zeros(Int64,nmstedges)
    tj = copy(ti)
    tv = zeros(Int64,nmstedges)
    k = 1
    for i = 1:nverts
        if pred[i] > 0
            j = pred[i]
            ti[k] = i
            tj[k] = j
            for rpi = rp[i]:rp[i+1]-1
                if ci[rpi] == j
                    tv[k] = ai[rpi]
                    break
                end
            end
            k = k + 1
        end
    end

    return (ti,tj,tv,nverts)
    
end

### Support more input:
mst_prim(A::MatrixNetwork) = mst_prim(A,false,1)
mst_prim(A::MatrixNetwork,full::Bool) = mst_prim(A,full,1)

### Support sparse matrices:

mst_prim(A::SparseMatrixCSC{T,Int64},full::Bool,u::Int64) where {T} = mst_prim(MatrixNetwork(A),full,u)
mst_prim(A::SparseMatrixCSC{T,Int64},full::Bool) where {T} = mst_prim(MatrixNetwork(A),full,1)
mst_prim(A::SparseMatrixCSC{T,Int64}) where {T} = mst_prim(MatrixNetwork(A),false,1)


### Output modifiers:
function mst_prim_matrix(A::MatrixNetwork,full::Bool,u::Int64)
    (ti,tj,tv,nverts) = mst_prim(A,full,u)
    T = sparse(ti,tj,tv,nverts,nverts)
    T = T + T'
    return T
end
mst_prim_matrix(A::MatrixNetwork) = mst_prim_matrix(A,false,1)
mst_prim_matrix(A::MatrixNetwork,full::Bool) = mst_prim_matrix(A,full,1)

mst_prim_matrix(A::SparseMatrixCSC{T,Int64},full::Bool,u::Int64) where {T} = mst_prim_matrix(MatrixNetwork(A),full,u)
mst_prim_matrix(A::SparseMatrixCSC{T,Int64},full::Bool) where {T} = mst_prim_matrix(MatrixNetwork(A),full,1)
mst_prim_matrix(A::SparseMatrixCSC{T,Int64}) where {T} = mst_prim_matrix(MatrixNetwork(A),false,1)
