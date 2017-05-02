"""
Ear Decomposition
-----------------
An ear decomposition of an undirected graph G is a partition of its set of edges into a
sequence of ears, such that the one or two endpoints of each ear belong to earlier ears
in the sequence and such that the internal vertices of each ear do not belong to any earlier ear.

The code is based on the algorithm given by Schmidt et.al in
"A simple test on 2-vertex- and 2-edge-connectivity" (Information Processing Letters Journal' 2013).

This approach is an important subroutine in many graph problems like Planarity Testing, Identifying the
Factor critical graphs, recognizing series parallel graphs, identifying Tri connected components etc.

"""

type Ear_decomposition_output
    ears::Array{Tuple{Int64,Int64,Int64},1} #disjoint ears or paths
    non_trivial_ears::Int64 #count of non trivial paths in the decomposition
    A::MatrixNetwork #MatrixNetwork
end

"""
`ear_decomposition!`
---
This function returns the number of non-trivial ears of the underlying
biconnected graph.

Functions
---------
- `ears = ear_decomposition!(A::MatrixNetwork,ears::Array{Tuple{Int64,Int64,Int64},1})`

Inputs
------
-`A`: The adjacency matrix.
`ears`: An empty array of tuples

Updates
-------
-`ears`: The ear decomposition of the graph.

Returns
-------
-`non_trivial_ears` : The count of non trivial ears in the adjacency matrix

Example
-------
A = load_matrix_network("ear_example")
B = MatrixNetwork(A)
ears = Array{Tuple{Int64,Int64,Int64},1}
non_trivial_ears = ear_decomposition!(A,ears)
"""

function ear_decomposition!(A::MatrixNetwork, ears::Array{Tuple{Int64,Int64,Int64},1})
    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    d=-1*ones(Int64,n)
    dt=-1*ones(Int64,n)
    rev_dt= -1*ones(Int64,n)
    ft=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    rs=Tuple{Int,Int}[]
    be = Array{Int64, 1}[]
    #be = Any[] #back edges array
    #ears=Tuple{Int64,Int64,Int64}[]
    ear_visited = zeros(Int64, A.n)
    for i=1:n
        push!(be, []) #Create an empty adjacencylist for back edges.
    end

    #start dfs
    t=1
    targethit=0
    for i=1:n
        v=i
        if d[v]>0
            continue
        end
        d[v]=1
        dt[v]=t
        rev_dt[t]=v;
	t=t+1
        ri=rp[v]
        push!(rs,(v,ri))
        while size(rs,1)>0
            (v,ri)=pop!(rs)
            while ri<rp[v+1]
                w=ci[ri]
                ri=ri+1
                if d[w]<0
                    d[w]=d[v]+1
                    pred[w]=v
                    push!(rs,(v,ri))
                    v=w
                    ri=rp[w]
                    dt[v]=t
		    rev_dt[t]=v
                    t=t+1
                elseif (pred[v]!=w)
		    if(dt[w]< dt[v])
		       push!(be[w], v) # Store the back edges encountered.
		    end
		end
            end
        end
        if full == 0
            break
        end
    end
    cn=1
    non_trivial=0
    for i=1:n
	startVertex = rev_dt[i]
	ear_visited[startVertex] = 1
	while !isempty(be[startVertex]) #A disjoint Path (or edge) is identified.
            flag=0
	    nextVertex = pop!(be[startVertex])
            push!(ears,(startVertex,nextVertex,cn))
            while (ear_visited[nextVertex] != 1) #Identify the ear(chain) and number it.
                flag=1
                ear_visited[nextVertex] = 1
                startVertex_n = nextVertex
                nextVertex = pred[startVertex_n]
                push!(ears,(startVertex_n,nextVertex,cn))
            end
            cn += 1
            if flag==1
                non_trivial+=1
            end
	end
    end
    return non_trivial
end

"""
`ear_decomposition`
---
This function returns the ear decomposition of the underlying biconnected graph.
It expects a biconnected graph as input

Functions
---------
- `ears = ear_decomposition(A::MatrixNetwork)

Inputs
------
-`A`: The adjacency matrix.

Returns
-------
-Returns a `Ear_decomposition_output` type that includes
`ears` : decomposition of biconnected component into a set of disjoint paths
`non_trivial_ears`: Count of non trivial ears in the adjacency matrix A

Example
-------
A = load_matrix_network("ear_example")
B = MatrixNetwork(A)
obj = ear_decomposition(B)
ears = obj.ears
count_non_trivial_ears = obj.non_trivial_ears
"""

function ear_decomposition(A::MatrixNetwork)
    ears=Tuple{Int64,Int64,Int64}[]
    non_trivial_ears=ear_decomposition!(A,ears)
    return (Ear_decomposition_output(ears,non_trivial_ears,A))
end


###############################
##    Conversion Functions    #
###############################

#CSC
ear_decomposition(A::SparseMatrixCSC) = biconnected_components(MatrixNetwork(A))

#Triplet
ear_decomposition(ei::Vector{Int64},ej::Vector{Int64}) = ear_decomposition(MatrixNetwork(ei,ej))

#CSR
ear_decomposition{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64) = ear_decomposition(MatrixNetwork(n,rp,ci,vals))
