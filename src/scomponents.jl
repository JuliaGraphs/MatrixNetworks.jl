
# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation
# require("MatrixNetworks.jl")
"""
Example
-------

- `cc` = scomponents(A)
- scomponents(A).num
- scomponents(A).sizes
- scomponents(A).map
- strong_components_map(A)     # if you just want the map
- enrich(scomponents(A)) # produce additional enriched output

Return information on the strongly connected components of a graph.
The method used in Tarjan's algorithm.
"""
:scomponents

###########################
##    Type Definitions    #
###########################
type strong_components_output
    map    # The indicator map
    sizes  # Array of sizes corresponding to map
    number # Int64
    A      # MatrixNetwork
end

type strong_components_rich_output
    reduction_matrix # the reduction matrix (restriction matrix)
    transitive_order
    transitive_map
end

####################
##    Functions    #
####################


"""
Return information on the strongly connected components of a graph that 
is the minimum required computation.
Example:
``MatrixNetworks.strong_components_map(MatrixNetworks.MatrixNetwork(sprand(5,4,0.5)))``	
"""
function strong_components_map(A::MatrixNetwork)

	# TODO, remove the dt variable here
	n=A.n
    sci=zeros(Int64,n)
    cn=1
    root=zeros(Int64,n)
    dt=zeros(Int64,n)
    t=0
    cs=zeros(Int64,n)
    css=0 # component stack
    rs=zeros(Int64,2*n)
    rp = A.rp
    ci = A.ci
    rss=0 # recursion stack holds two nums (v,ri)
    
    # start dfs at 1
    for sv=1:n
        v=sv
        if root[v]>0
            continue
        end
        rss=rss+1
        rs[2*rss-1]=v
        rs[2*rss]=rp[v] #add v to the stack
        root[v]=v 
        sci[v]=-1
        dt[v]=t
        t=t+1
        css=css+1
        cs[css]=v #add w to component stack
        while rss>0
            v=rs[2*rss-1] # pop v from the stack
            ri=rs[2*rss]
            rss=rss-1 
            while ri<rp[v+1]
                w=ci[ri]
                ri=ri+1
                if root[w]==0
                    root[w]=w; sci[w]=-1;  dt[w]=t; t=t+1;
                    css=css+1 # add w to component stack
                    cs[css]=w 
                    rss=rss+1
                    rs[2*rss-1]=v
                    rs[2*rss]=ri # add v to the stack
                    v=w
                    ri=rp[w]
                    continue
                end
            end
            for ri=rp[v]:rp[v+1]-1
                w=ci[ri]
                if sci[w]==-1
                    if dt[root[v]]>dt[root[w]]
                       root[v]=root[w]
                    end
                end
            end
            if root[v]==v
                while css>0
                    w=cs[css]
                    css=css-1 
                    sci[w]=cn
                    if w==v
                        break
                    end
                end
                cn=cn+1
            end
        end
    end
    return sci
end

###############################
##    Conversion Functions    #
###############################

"""
Example
A = sprand(5,5,0.5)
MatrixNetworks.strong_components_map(A)
Example:
ei = [1;2;3]
ej = [2;4;1]
MatrixNetworks.strong_components_map(ei,ej)
"""

strong_components_map(A::SparseMatrixCSC{Float64,Int64}) = scomponents(MatrixNetwork(A))
strong_components_map(ei,ej) = strong_components_map(MatrixNetwork(ei,ej))



function scomponents(A::MatrixNetwork)
    map = strong_components_map(A)
	number = maximum(map)
	sizes = zeros(Int64,number);
    for i = 1:length(map)
        sizes[map[i]] += 1;
    end
    
    return strong_components_output(map, sizes, number, A)
end

###############################
##    Conversion Functions    #
###############################
# TODO: double check output of enrich and scomponents
# scomponents(A::SparseMatrixCSC{Float64,Int64}) = scomponents(MatrixNetwork(A))
# scomponents(ei,ej) = scomponents(MatrixNetwork(ei,ej))


""" 
This function adds the following helpers variables
* reduction_matrix - a reduction matrix to project down to the component graph
* transitive_order - a transitive ordering of the components
* transitive_map - a map to components that respects the transitive ordering
* largest
"""
:enrich
function enrich(rval::strong_components_output)
    ci = rval.map
    sizes = rval.sizes
    ncomp = maximum(ci)
    R = sparse([1:size(A,1);],ci,1,size(A,1),ncomp)
    CG = R'*A*R
    return strong_components_rich_output(R,CG,CG)
end





