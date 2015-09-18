#TODO: documentation
#TODO: provide example and documentation

"""
DIJKSTRA documentation here

Example
-------

# Find the minimum travel time between Los Angeles (LAX) and
# Rochester Minnesota (RST).

(A,xy,labels) = load_matrix_network_metadata("airports");\n
A = -A; # fix funny encoding of airport data\n
lax = 247; rst = 355;\n
(d,pred) = dijkstra(A,lax);\n
@printf("Minimum time: %d\n",d[rst]); #Print the path\n
@printf("Path:\n");\n
u = rst;\n
while(u != lax)\n
    @printf("%s <-- ", labels[u])\n
    u = pred[u];\n
    if (u == lax)\n
        @printf("%s\n", labels[lax])\n
    end\n
end
"""

function dijkstra(A::MatrixNetwork,u::Int64)
    (rp,ci,ai) = (A.rp, A.ci, A.vals)
    return dijkstra_internal(rp,ci,ai,u)
end

function dijkstra{F}(A::SparseMatrixCSC{F,Int64},u::Int64)
    (rp,ci,ai) = sparse_to_csr(A)
    return dijkstra_internal(rp,ci,ai,u)
end

function dijkstra_internal{F}(rp::Vector{Int64},ci::Vector{Int64},ai::Vector{F},u::Int64)

    if any(ai.<0)
        error("dijkstra''s algorithm cannot handle negative edge weights")
    end
    
    n = length(rp) - 1
    d = Inf*ones(Float64,n)
    T = zeros(Int64,n)
    L = zeros(Int64,n)
    pred = zeros(Int64,length(rp)-1)
    
    n = 1
    T[n] = u
    L[u] = n # oops, n is now the size of the heap
    
    # enter the main dijkstra loop
    d[u] = 0
    while n > 0
        v = T[1]
        ntop = T[n]
        T[1] = ntop
        L[ntop] = 1
        n = n - 1 # pop the head off the heap
        k = 1
        kt = ntop # move element T[1] down the heap
        while true
            i = 2*k
            if i > n
                break
            end      # end of heap
            if i == n
                it = T[i]        # only one child, so skip
            else               # pick the smallest child
                lc = T[i]
                rc = T[i+1]
                it = lc
                if d[rc] < d[lc]
                    i = i+1
                    it = rc
                end # right child is smaller
            end
            if d[kt] < d[it]
                break   # at correct place, so end
            else 
                T[k] = it
                L[it] = k
                T[i] = kt
                L[kt] = i
                k=i # swap
            end
        end       # end heap down
        # for each vertex adjacent to v, relax it
        for ei = rp[v]:rp[v+1]-1       # ei is the edge index
            w = ci[ei]
            ew = ai[ei]          # w is the target, ew is the edge weight
            # relax edge (v,w,ew)
            if d[w] > d[v] + ew
                d[w] = d[v] + ew
                pred[w] = v
                # check if w is in the heap
                k = L[w]
                onlyup=false
                if k == 0
                    # element not in heap, only move the element up the heap
                    n = n+1
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
                    i = 2 * k
                    if i > n
                        break
                    end          # end of heap
                    if i == n
                        it = T[i]    # only one child, so skip
                    else            # pick the smallest child
                        lc = T[i]
                        rc = T[i+1]
                        it = lc
                        if d[rc] < d[lc]
                            i = i+1
                            it = rc
                        end # right child is smaller
                    end
                    if d[kt] < d[it]
                        break      # at correct place, so end
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
                    j2 = round(Int,floor(j/2))
                    tj2 = T[j2]    # parent element
                    if d[tj2] < d[tj]
                        break      # parent is smaller, so done
                    else                         # parent is larger, so swap
                        T[j2] = tj
                        L[tj] = j2
                        T[j] = tj2
                        L[tj2] = j
                        j=j2
                    end
                end
            end
        end
    end
    return (d,pred)
end

