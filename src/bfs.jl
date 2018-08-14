"""
BFS
---
    compute breadth first search distances and returns a distance(d), 
    the discover time(dt), predecessor array(pred) in the tuple (d,dt,pred)
    pred[i] = 0 if vertex i is in a component not reachable from u and i != u.
    Search stops when it hits the vertex target.

Functions
---------
* (d,dt,pred) = bfs(A::MatrixNetwork,u::Int64,target::Int64)
* (d,dt,pred) = bfs{T}(A::SparseMatrixCSC{T,Int64}),u::Int64,target::Int64)
* (d,dt,pred) = bfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64,target::Int64)\n
If target is not specified, it is assigned to 0

Example
-------
~~~
A = load_matrix_network("bfs_example")
(d,dt,pred) = bfs(A,1)
~~~
"""
function bfs(A::MatrixNetwork,u::Int64,target::Int64)
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
                if w==target
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
function bfs(A::SparseMatrixCSC{T,Int64},u::Int64) where T
    return bfs(MatrixNetwork(A),u)
end

function bfs(A::SparseMatrixCSC{T,Int64},u::Int64,target::Int64) where T
    return bfs(MatrixNetwork(A),u,target)
end

## Triplet Format:
function bfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64)
    return bfs(MatrixNetwork(ei,ej),u)
end

function bfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64,target::Int64)
    return bfs(MatrixNetwork(ei,ej),u,target)
end

## CSR sparse matrices:
function bfs(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64,u::Int64) where T
    return bfs(MatrixNetwork(n,rp,ci,vals),u)
end

function bfs(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64,u::Int64,target::Int64) where T
    return bfs(MatrixNetwork(n,rp,ci,vals),u,target)
end