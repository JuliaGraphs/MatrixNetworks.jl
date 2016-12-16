"""

biconnected Component Decomposition
----------------------------------

Any connected graph decomposes into a tree of biconnected components called the block-cut tree of the graph. 
The blocks are attached to each other at shared vertices called cut vertices or articulation points.

This implementation is based on the algorithm provided by Tarjan in "Depth-First Search and Linear Graph Algorithms".  


'1' in the entry of the Articulation point array denotes that the vertex is an articulation point.


Functions
---------
"""

type biconnected_components_output
    bcc_edges::Array{Tuple,1} #list of bcc edges
    Articulation_points::Vector{Int64}  
    A::MatrixNetwork #MatrixNetwork
end



function biconnected_component(A::MatrixNetwork)
    n=A.n
    rp=A.rp
    ci=A.ci
    bcc_edges=Tuple{Int,Int,Int}[]
    cn=1
    low=zeros(Int64,n)
    dt=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    t=0
    largest_bcc=0
    edge_count_bcc=0
    rs=Tuple{Int,Int,Int}[]
    cs=Tuple{Int,Int}[]
    root_children= 0
    articulation = zeros(Int64,n)
    
    #start dfs at 1
    for sv=1:n
        v=sv
        if dt[v]>0
            continue
        end
        root_children = 0
        low[v]=t
        dt[v]=t;
        t=t+1;
        w = ci[rp[v]]
        if (v == w) && (rp[v]+1 == rp[v+1])
            push!(bcc_edges, (v,v,cn))
            cn = cn+1
            articulation[v]=1
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
                        push!(cs, (parent,child))
                    end
                else
                    low[child] = dt[child] = t
                    t+=1
                    push!(rs, (parent, child, rp[child]))
                    push!(cs, (parent,child))
                end
            else
                pop!(rs)
                if size(rs,1) > 1
                    if low[parent] >= dt[g_parent]  # g_parent is an Articulation point
                        articulation[g_parent]=1
                        while true
                            (u,v)=pop!(cs)
                            push!(bcc_edges, (u,v,cn))  #Assign component number to the corresponding edges
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
                        (u,v) = pop!(cs)
                        push!(bcc_edges,(u,v,cn))
                        if u == g_parent && v == parent
                            break
                        end
                    end
                    cn+=1
                end
            end
        end
    end	
    return (bcc_edges, articulation, cn)
end
 


function biconnected(A::MatrixNetwork)
    (mapping,Articulation_points) = biconnected_component(A)
    return  biconnected_components_output(mapping,Articulation_points,A)
end


function rich_output(A::MatrixNetwork, biconnected_component_id::Int64 = 0)
    (mapping,Articulation_points,cn) = biconnected_component(A)
    component_edges = Tuple{Int,Int}[]
    max_id=max_count=current_count=0
    prev_id=0
    if biconnected_component_id <= 0
        for (u,v,id) in mapping
            current_id=id
            if current_id == prev_id
                current_count+=1
            else
                if max_count < current_count
                    max_count = current_count
                    max_id = prev_id
                end
                current_count = 1
                prev_id = current_id
            end      
        end
        biconnected_component_id = max_id
    end
    
    if biconnected_component_id >= cn 
        return component_edges
    end           
    for (u,v,id) in mapping
        if id == biconnected_component_id
            push!(component_edges,(u,v))
        end
    end
    return component_edges
end
