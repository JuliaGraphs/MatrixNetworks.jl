"""
biconnected component decomposition
----------------------------------

Any connected graph decomposes into a tree of biconnected components 
called the block-cut tree of the graph. The blocks are attached to each 
other at shared vertices called cut vertices or articulation points.

This implementation is based on the algorithm provided by Tarjan 
in "Depth-First Search and Linear Graph Algorithms".  
"""
mutable struct Biconnected_components_output
    map::Vector{Int64} #biconnected_component_number
    articulation_points::Vector{Bool}
    number::Int64
    A::MatrixNetwork #MatrixNetwork
end

"""
`biconnected_components!`
---
This function returns the number of biconnected components in the 
underlying graph. It expects an undirected graph as its input.

Functions
---------
- `number = biconnected_components!(A::MatrixNetwork, articulation::Vector{Bool}, map::Vector{Int64})`

Inputs
------
- `A`: the adjacency matrix.
- `articulation`: A boolean array, where each element is initialized to false.
- `map`: Vector of size equal to the number of edges.

Returns
-------
- `cn`: The number of biconnected components in the graph

Example
-------
A = load_matrix_network("biconnected_example")
B = MatrixNetwork(A)
number_of_components = biconnected_components!(B, zeros(Bool,0), zeros(Int64,0))
"""
function biconnected_components!(A::MatrixNetwork, articulation::Vector{Bool}, map::Vector{Int64})
    n=length(A.rp)-1
    rp=A.rp
    ci=A.ci
    components = length(map) >= length(ci) 
    cn=1
    low=zeros(Int64,n)
    dt=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    t=0
    largest_bcc=0
    edge_count_bcc=0
    rs=Tuple{Int,Int,Int}[]
    cs=Tuple{Int,Int,Int}[]
    root_children= 0
    art = length(articulation) >= n 

    #start dfs at 1.
    for sv=1:n
        v=sv
        if dt[v]>0
            continue
        end
        root_children = 0
        low[v]=t
        dt[v]=t
        t=t+1
        w = ci[rp[v]]
        if (v == w) && (rp[v]+1 == rp[v+1])
            if components 
                map[rp[v]]=cn
            end
            cn = cn+1
            if art 
                articulation[v]=1 
            end
        else
            push!(rs, (v,v,rp[v])) #add v,rp[v] to stack
        end

        while size(rs, 1) > 0
            (g_parent, parent, children_index) = rs[end]
            child_index = children_index
            if child_index != rp[parent+1]
               child = ci[child_index]
               children_index = children_index+1
               (g_parent, parent, temp) = pop!(rs)
               push!(rs, (g_parent, parent, children_index))
               if (g_parent == child) || (parent == child)		
                   continue   # skip tree edges and self loops
               end
               if dt[child] >= 0
                    if dt[child] <= dt[parent]
                        if dt[child] < low[parent] #Update Low
                            low[parent] = dt[child]
                        end
                        push!(cs, (parent,child,child_index))
                    end
                else
                    low[child] = dt[child] = t
                    t+=1
                    push!(rs, (parent, child, rp[child]))
                    push!(cs, (parent,child,child_index))
                end
            else
                pop!(rs)
                if size(rs,1) > 1
                    if low[parent] >= dt[g_parent]  # g_parent is an Articulation point
                        if art 
                            articulation[g_parent]=1
                        end
                        while true
                            (u,v,child_index)=pop!(cs)
                            if components
                                map[child_index]=cn #assign the component number of the corresponding edge
                            end
                            if u == g_parent && v == parent
                               break
                            end
                        end
                        cn+=1
                    end
                    if low[parent] < low[g_parent]
                       low[g_parent] = low[parent]
                    end
                elseif size(rs,1) == 1
                    root_children +=1
                    while true
                        (u,v,child_index) = pop!(cs)
                        if components
                            map[child_index]=cn
                        end
                        if u == g_parent && v == parent
                            break
                        end
                    end
                    cn+=1
                end
            end
        end
    end
    cn=cn-1
    
    return cn
end


"""
`biconnected_components`
-----------------------
This function requires a symmetric matrix as input. Depending on the user input this
function returns either the biconnected component number associated with each edge 
or articlation points or both.

Inputs
---------
- `A`: the adjacency matrix
- Optional Keyword inputs
  - `art=true`: returns the articulation points of the graph.
  - `components=true`:returns the biconnected component labels associated with each 
  edge.

Returns
-------
- Returns a `Biconnected_components_output` type which includes
`map` : biconnected component labels associated with each edge, 
`articulation_points`: boolean array that signifies whether a vertex is an articulation point and
`number`: Number of biconnected components in the graph.

Example
-------
A = load_matrix_network("biconnected_example")
B = MatrixNetwork(A)
bcc = biconnected_components(B)
map = bcc.map
articulation_vector = bcc.articulation_points
number_of_components = bcc.number
"""
function biconnected_components(A::MatrixNetwork; art::Bool = true, components::Bool = true)
    map = components ? zeros(Int64, length(A.ci)) : zeros(Int64, 0)
    articulation = art ? zeros(Bool, A.n) : zeros(Bool, 0)
    cn = biconnected_components!(A, articulation, map)
    return Biconnected_components_output(map,articulation,cn,A)
end

###############################
##    Conversion Functions    #
###############################

#CSC
biconnected_components(A::SparseMatrixCSC;kwargs...) = biconnected_components(MatrixNetwork(A);kwargs...)

#Triplet
biconnected_components(ei::Vector{Int64},ej::Vector{Int64};kwargs...) = biconnected_components(MatrixNetwork(ei,ej);kwargs...)

#CSR
biconnected_components(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64; kwargs...) where {T} = biconnected_components(MatrixNetwork(n,rp,ci,vals);kwargs...)
