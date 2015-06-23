# TODO: more testing and check documentation
# TODO: add more examples
include("MatrixNetwork.jl")
"""
BFS compute breadth first search distances and returns a distance (d), 
the discover time(dt), predecessor array (pred) in the tuple (d,dt,pred).
pred[i] = 0 if vertex i is in a component not reachable from u and i != u.
Example:\n
(d,dt,pred) = bfs(A,u,v) # search stops when it hits the vertex v\n
(d,dt,pred) = bfs(A,u)\n

Example
-------
file_path = Pkg.dir("MatrixNetworks/data/bfs_example.smat")\n
A = readSMAT(file_path)\n
(d,dt,pred) = bfs(MatrixNetwork(A),1)\n
OR: (d,dt,pred) = bfs(A,1)
"""

function bfs(A::MatrixNetwork,u::Int64,v::Int64)
    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    d=-1*ones(Int64,n)
    dt=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    sq=zeros(Int64,n)
    sqt=0
    sqh=0 # search queue and search queue tail/head
    
    # start bfs at u
    sqt=sqt+1
    sq[sqt]=u
    t=0
    d[u]=0
    dt[u]=t
    t=t+1
    pred[u]=u
    while sqt-sqh>0
        sqh=sqh+1
        v=sq[sqh] # pop v off the head of the queue
        for ri=rp[v]:rp[v+1]-1
            w=ci[ri]
            if d[w]<0
                sqt=sqt+1
                sq[sqt]=w
                d[w]=d[v]+1
                dt[w]=t
                t=t+1
                pred[w]=v
                if w==v
                    return (d,dt,pred)
                end
            end
        end
    end
    return (d,dt,pred)
end

############################
### Additional functions ###
############################

function bfs(A::MatrixNetwork,u::Int64)
    return bfs(A,u,0)
end

## CSC sparse matrices:
function bfs(A::SparseMatrixCSC{Float64,Int64})
    retrun bfs(MatrixNetwork(A))
end

function bfs(A::SparseMatrixCSC{Float64,Int64},u::Int64)
    retrun bfs(MatrixNetwork(A),u)
end

## Triplet Format:
function bfs(ei::Vector{Int64},ej::Vector{Int64})
    retrun bfs(MatrixNetwork(ei,ej))
end

function bfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64)
    retrun bfs(MatrixNetwork(ei,ej),u)
end

## CSR sparse matrices:
function bfs(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{Float64},n::Int64)
    retrun bfs(MatrixNetwork(n,rp,ci,vals))
end

function bfs(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{Float64},n::Int64,u::Int64)
    retrun bfs(MatrixNetwork(n,rp,ci,vals),u)
end