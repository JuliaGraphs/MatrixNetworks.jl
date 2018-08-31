"""
FLOYDWARSHALL
-------------
    compute all shortest paths using the Floyd-Warshall algorithm.
    
    (D,P) = floydwarshall(A) returns the shortest distance matrix between all pairs
    of nodes in the graph A in matrix D.  If A has a negative weight cycle, then this
    algorithm will throw an error. P is the matrix of predecessors.

Functions
---------
- (D,P) = floydwarshall(A::MatrixNetwork)
- (D,P) = floydwarshall{T}(A::SparseMatrixCSC{T,Int64})

Example
-------
~~~
A = load_matrix_network("all_shortest_paths_example")
(D,P) = floydwarshall(A)
~~~
"""
:floydwarshall

## setup functions:

function floydwarshall_phase1(A::MatrixNetwork)
    (nzi,nzj,nzv) = csr_to_sparse(A.rp,A.ci,A.vals)
    return (nzi,nzj,nzv,A.n)
end

function floydwarshall_phase1(A::SparseMatrixCSC{T,Int64}) where T
    (ri,ci,ai) = findnz(A)
    return (ri,ci,ai,A.n)
end

function floydwarshall_phase2(ri::Vector{Int64},ci::Vector{Int64},ai::Vector{T},n::Int64) where T

    nz = length(ai)
    D = Inf*ones(Int64,n,n)
    
    #TODO: check: always compute P or give the option of just computing D?
    
    P = zeros(Int64,n,n)
    # initialize the distance and predecessor matrix
    for ei = 1:nz
        i = ri[ei]
        j = ci[ei]
        v = ai[ei]
        if v < D[i,j]
            D[i,j] = v
            P[i,j] = i
        end
    end
    
    ids = (LinearIndices((n,n)))[CartesianIndex.(1:n, 1:n)]
    D[ids] .= 0 # set diagonal to 0
    
    for k=1:n
        for i=1:n
            for j=1:n
                if D[i,k]+D[k,j] < D[i,j]
                    D[i,j] = D[i,k]+D[k,j]
                    P[i,j] = P[k,j]
                end
            end
        end
    end
    
    if any(diag(D).<0)
        warn("floydwarshall:negativeCycle","negative weight cycle detected")
    end
    
    return (D,P)
end


## floyd warshall
function floydwarshall(A::MatrixNetwork)
    (nzi,nzj,nzv,n) = floydwarshall_phase1(A)
    (D,P) = floydwarshall_phase2(nzi,nzj,nzv,n)
    return (D,P)
end

function floydwarshall(A::SparseMatrixCSC{T,Int64}) where T
    (nzi,nzj,nzv,n) = floydwarshall_phase1(A)
    (D,P) = floydwarshall_phase2(nzi,nzj,nzv,n)
    return (D,P)
end
