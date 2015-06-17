# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation
# TODO: add more examples
# TODO support struct

"""
Example
-------

using MAT

file_path = Pkg.dir("MatrixNetworks/data/clique-10.mat")
  
file = matopen(file_path)

A = read(file,"A")

close(file)

cc = clustercoeffs(MatrixNetwork(A))

clustercoeffs compute undirected clustering coefficients for a graph
clustercoeffs(A) computes a normalized, weighted clustering coefficients from a graph
represented by a symmetric adjacency matrix A.
clustercoeffs(A,weighted,normalized), with weighted and normalized boolean values indicate
whether the computation has to be weighted and/or normalized.
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
	# TODO: fix this condition: Maybe add one more field to MatrixNetwork to save computation  
	# because we're already performing a transpose in the very beginning  
	#     if !(A.a == A'.a)
	#         #TODO: a more descriptive error stmt check if I can refer to another fn
	#         error("Only undirected (symmetric) inputs are allowed")
	#     end
    
    # if usew
    #     (rp,ci,ai)=sparse_to_csr(A)
    # else
    #     (rp,ci) = sparse_to_csr(A)
    # end
    # in julia's implementation of sparse to csr, it returns the same output so:
    #(rp,ci,ai)=sparse_to_csr(A)
    # passing a MatrixNetwork instead
    
    (rp,ci,ai) = (A.rp,A.ci,A.vals);
    if length(find(ai.<0)) != 0
        # TODO: fix error print statement
        error("only positive edge weights allowed")
    end
    
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