# TODO: more testing and better documentation

"""
FLOYDWARSHALL Compute all shortest paths using the Floyd-Warshall algorithm.

(D,P) = floydwarshall(A) returns the shortest distance matrix between all pairs
of nodes in the graph A in matrix D.  If A has a negative weight cycle, then this
algorithm will throw an error. P is the matrix of predecessors.

Example
-------
file_path = Pkg.dir("MatrixNetworks/data/all_shortest_paths_example.smat")\n
A = readSMAT(file_path)\n
(D,P) = floydwarshall(MatrixNetwork(A))\n
"""

function floydwarshall(A::MatrixNetwork)

    (rp,ci,ai) = (A.rp,A.ci,A.vals)
    nz = length(ai)
    n = A.n
    D = Inf*ones(Float64,n,n)
    
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
    
    ids = sub2ind(n,1:n,1:n)
    D[ids] = 0 # set diagonal to 0

#     for i=1:n
#         D[i,i] = 0
#     end # set diagonal to 0
    
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
    
    if any(diag(D))<0
        warn("floydwarshall:negativeCycle","negative weight cycle detected")
    end
    
    return (D,P)
end
