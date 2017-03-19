"""

biconnected Component Decomposition
----------------------------------

Any connected graph decomposes into a tree of biconnected components 
called the block-cut tree of the graph. The blocks are attached to each 
other at shared vertices called cut vertices or articulation points.

This implementation is based on the algorithm provided by Tarjan 
in "Depth-First Search and Linear Graph Algorithms".  


Functions
---------

biconnected_components(A::MatrixNetwork, art::Bool, components::Bool)
--------------------------------------------------------------------
Returns an array that represents the biconnected component number for the corrsponding 
edge in its csr representation. Set the art or components parameters to true 
depending on whether the user requires only articulation points or biconnected components 
or both. '1' in the entry of the articulation point array denotes that the vertex is an articulation point.

enrich(A::MatrixNetwork, biconnected_component_id::Int64)
---------------------------------------------------------------------
Returns the ordered list of edges part of each biconnected component.


Example
-------

A = load_matrix_network("biconnected_example")
B = MatrixNetwork(A)
bcc = biconnected(B)
bcc.map
bcc.articulation_points

"""

type biconnected_components_output
    map::Vector{Int64} #biconnected_component_number
    articulation_points::Vector{Int64}  
    components::Int64
    A::MatrixNetwork #MatrixNetwork
end

type biconnected_components_rich_output
    map::Vector{Tuple{Int,Int,Int}} #ordered list of biconnected_component_numbers
end

function biconnected_components(A::MatrixNetwork, art::Bool = true, components::Bool = true)
    n=length(A.rp)-1
    rp=A.rp
    ci=A.ci
    map=components ? zeros(Int64,length(ci)) : zeros(Int64,0)
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
    articulation=art ? zeros(Int64,n) : zeros(Int64,0)

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
    return biconnected_components_output(map,articulation,cn,A)
end

function biconnected_components!(obj::biconnected_components_output)
    return obj.components
end

function enrich_helper(A::MatrixNetwork, mapping::Vector{Int64})  #Initializer function for enrich
    bcc_edges = Tuple{Int,Int,Int}[]
    n = length(A.rp)-1
    ci = A.ci
    rp = A.rp
    for i = 1:n
        for j = rp[i]:rp[i+1]-1
            if mapping[j]!=0
                push!(bcc_edges, (i,ci[j],mapping[j]))
            end
        end 
    end
    bcc_edges = sort(bcc_edges, by = x -> x[3])
    return bcc_edges
end

function enrich(obj::biconnected_components_output)
    bcc_edges = enrich_helper(obj.A,obj.map)    #Ordered list of biconnected component edges
    return biconnected_components_rich_output(bcc_edges)
end


###############################
##    Conversion Functions    #
###############################

#CSC
biconnected_components{T}(A::SparseMatrixCSC{T,Int64}, art::Bool = true, component::Bool = true) = biconnected_components(MatrixNetwork(A), art, component)
#Triplet
biconnected_components(ei::Vector{Int64},ej::Vector{Int64},art::Bool = true, component::Bool = true) = biconnected_components(MatrixNetwork(ei,ej), art, component)
#CSR
biconnected_components{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64, art::Bool = true, component::Bool = true) = biconnected_components(MatrixNetwork(n,rp,ci,vals), art, component)

#CSC
enrich{T}(A::SparseMatrixCSC{T,Int64}, biconnected_component_id::Int64 = 0) = enrich(MatrixNetwork(A), biconnected_component_id)
#Triplet
enrich(ei::Vector{Int64},ej::Vector{Int64}, biconnected_component_id::Int64 = 0) = enrich(MatrixNetwork(ei,ej), biconnected_component_id)
#CSR
enrich{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64, biconnected_component_id::Int64 = 0) = enrich(MatrixNetwork(n,rp,ci,vals), biconnected_component_id)
