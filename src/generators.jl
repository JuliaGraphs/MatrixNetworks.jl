
"""
`erdos_renyi_undirected`
========================

Generate an undirected Erdős-Rényi graph.

**The current implementation uses sprand, this may change in the future.**
**Do not depend on this routine for reliable output between versions.**


Input
-----
- `n`: the number of nodes
- `p`: the probability of an edge, or the average degree. 
   if \$`p` >= 1\$, then \$`p`\$ is interpreted as an average degree
   instead of a probability. (There is no point in generating
   an Erdős-Rényi graph with probability \$`p`=1\$)

Output
------
- A matrix network type for the Erdős-Rényi graph.

Example
-------
~~~~ 
# show the connected phase transition
n = 100
avgdegs = linspace(1.,2*log(n),100) 
compsizes = map( (dbar) -> 
        maximum(scomponents(erdos_renyi_undirected(n,dbar)).sizes),
    avgdegs )
using Plots
unicodeplots()
plot(avgdegs,compsizes,xaxis=("average degree"),yaxis=("largest component size"))    
~~~~    
"""
function erdos_renyi_undirected(n::Int, p::Float64)
    if p < 0 || p > n throw(DomainError()) end
    
    if p >= 1. # interpret as average degree
        p = p/n # convert to probability
    end
    A = sprand(n,n,p)
    Aup = triu(A,1)
    Asym = max(Aup,Aup')
    return _matrix_network_direct(Asym,1)
end

erdős_rényi_undirected(n,p) = erdos_renyi_undirected(n,p) 

"""
`erdos_renyi_directed`
========================

Generate an directed Erdős-Rényi graph.

**The current implementation uses sprand, this may change in the future.**
**Do not depend on this routine for reliable output between versions.**


Input
-----
- `n`: the number of nodes
- `p`: the probability of an edge. 

Output
------
- A matrix network type for the Erdős-Rényi graph. 
"""
function erdos_renyi_directed(n::Int, p::Float64)
    if p < 0 || p > n throw(DomainError()) end
    if p >= 1. # interpret as average degree
        p = p/n # convert to probability
    end
    
    A = sprand(n,n,p)
    
    return _matrix_network_direct(A-diag(diag(A)),1) # directions don't matter
end

erdős_rényi_directed(n,p) = erdos_renyi_directed(n,p)

"""
`chung_lu_undirected`
========================

Generate an approximate undirected Chung-Lu graph. The approximation is because we draw
exactly |E| edges where each edge is sampled from the Chung-Lu model. But then
we discard duplicate edges and self-loops. So the new graph will always have fewer
edges than the input degree sequence.

**This will likely change in future versions and provide an exact Chung-Lu model.**

If the graph is 

Usage
-----
- `chung_lu_undirected(d)`
- `chung_lu_undirected(d,nedges)` 

Input
-----
- `d`: the degree sequence vector. This vector is non-negative and has the expected
degree for each vertex.

Output
------
- A MatrixNetwork for the undirected graph that results from the Chung-Lu sample.

Example
-------
~~~~
A = load_matrix_network("tapir")
d = vec(sum(A,1))
B = sparse(chung_lu_undirected(d))
nnz(A)
nnz(B)
~~~~

"""
function chung_lu_undirected end

function chung_lu_undirected(d::Vector{Int})
    chung_lu_undirected(d, floor(Int,sum(d)/2))
end

function chung_lu_undirected(d::Vector{Int}, nedges::Int)
    n = length(d)

    # TODO find some standard function for this
    for v in d
        if v < 0
            throw(DomainError())
        end
    end
    
    if nedges < 0 || nedges > div(n*(n-1),2)
        throw(ArgumentError("nedges $nedges is too large for $n node undirected graph")) 
    end
    
    nodevec = zeros(Int,sum(d))
    curedge = 1
    for i=1:n
        for j=1:d[i]
            nodevec[curedge] = i
            curedge += 1
        end
    end
    
    ei, ej = unique_edge_sample_undirected(nodevec, nedges)
    A = sparse(ei,ej,1.,n,n)

    return _matrix_network_direct(A) # avoid the transpose  
end

function unique_edge_sample_undirected(nodevec, nedges::Int)
    edges = Set{Tuple{Int,Int}}()
    sizehint!(edges, nedges)
        
    curedges = 0
    while curedges < nedges
        src = rand(nodevec)
        dst = rand(nodevec)
        if src == dst
            continue
        else
            # fix the order 
            if src < dst
                dst, src = src, dst
            end
            if !((src,dst) in edges)
                push!(edges, (src,dst))
                curedges += 1
            end
        end
    end

    ei = zeros(Int,nedges*2)
    ej = zeros(Int,nedges*2)
    
    k = 1
    for edge in edges
        ei[k] = edge[1]
        ej[k] = edge[2]
        k+=1    
        ei[k] = edge[2]
        ej[k] = edge[1]
        k+=1
    end
    
    return ei, ej
end

"""
Not public right now
"""
function _chung_lu_dense_undirected(d::Vector{Int})
    M = zeros(Int,n,n)
    idenom = 1/sum(d)
    for i=1:n
        for j=i+1:n
            if rand(Float64) < d[i]*d[j]*idenom
                M[i,j] = 1
            end
        end
    end
    M = M + M'
    return _matrix_network_direct(sparse(M))
end

"""
The internal Havel-Hakimi function has an optional store
behavior that saves the edges as they come out of the algorithm.
This enables us to generate a Havel Hakimi graph, which can
be useful.
"""
function _havel_hakimi(degs::Vector{Int}, store::Bool, ei::Vector{Int}, ej::Vector{Int})
    q = Collections.PriorityQueue(Int,Int,Base.Order.Reverse)
    n = length(degs)
    effective_n = n
    degsum = 0
    dmax = 0
    for (i,d) in enumerate(degs)
        q[i] = d
        degsum += d
        dmax = max(d,dmax)
        if d < 0
            throw(ArgumentError("the degree sequence must be non-negative"))
        end
    end
    
    if mod(degsum,2) != 0; return false; end
    if n > 0 && dmax >= n; return false; end # n > 0 checks for the empty graph
    
    if store
        resize!(ei,degsum)
        resize!(ej,degsum)
    end
    
    dlist = Vector{Pair{Int,Int}}(dmax)
    enum = 1
    
    while !isempty(q)
        vi,d = Collections.peek(q) # vi is the cur vertex, d is the cur deg
        Collections.dequeue!(q)    # remove it
        for n=1:d                  # make a list of each neighbor
            if isempty(q); return false; end
            dlist[n] = Collections.peek(q)
            Collections.dequeue!(q)
        end
        # now "add" an edge from vi->neighbor, and thus, decrease it's 
        # degree when we re-add it
        for n=1:d
            if store
                ei[enum] = vi
                ej[enum] = dlist[n][1]
                enum += 1
                ei[enum] = dlist[n][1]
                ej[enum] = vi
                enum += 1
            end
            if dlist[n][2] < 1
                return false
            elseif dlist[n][2] == 1
                # don't both re-adding the vertex
            else
                q[dlist[n][1]] = dlist[n][2] - 1
            end
        end
    end
    return true
end

"""
`is_graphical_sequence`
=======================

Check whether or not a degree sequence is graphical,
which means that it is a valid degree sequence for 
an undirected graph.

Note that this does not mean it is a valid degree
sequence for a connected undirected graph. So,
for instance, 
`[1,1,1,1]` 
is a valid degree sequence for two
disconnected edges

Usage
-----
`is_graphical_sequence(d)` returns true or false 

Input
-----
- `d::Vector{Int}`:  a vector of integer valued degrees

Output
------
- a boolean that is true if the sequence is graphical
""" 
function is_graphical_sequence(d::Vector{Int})
    return _havel_hakimi(d, false, Int[],Int[])
end


"""
`havel_hakimi_graph`
====================

Create a graph with a given degree sequence 

Usage
-----
`A = havel_hakimi_graph(d)` returns an instance of the 
a graph with degree sequence d or throws ArgumentError
if the degree sequence is not graphical.   

Input
-----
- `d::Vector{Int}`:  a vector of integer valued degrees

Output
------
-`A`: a matrix network for the undirected graph that
results from the Havel-Hakimi procedure.
""" 
function havel_hakimi_graph(d::Vector{Int})
    ei = Int[]
    ej = Int[]
    if _havel_hakimi(d, true, ei, ej) == false
        throw(ArgumentError("the degree sequence is not graphical"))
    end
    return MatrixNetwork(ei,ej,length(d))
end

# TODO Add chung-lu for general floating point weights
# via union of ER graphs add 