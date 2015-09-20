"""
MST_PRIM
--------
compute a minimum spanning tree with Prim's algorithm.\n
T = mst_prim_matrix(A) computes a minimum spanning tree T using Prim's algorithm
for the spanning tree of a graph with non-negative edge weights.

T = mst_prim_matrix(A,false) produces an MST for just the component at A containing
vertex 1.  T = mst_prim_matrix(A,0,u) produces the MST for the component
containing vertex u.

(ti,tj,tv,nverts) = mst_prim(A) returns the edges from the matrix and does not
convert to a sparse matrix structure.  This saves a bit of work and is
required when there are 0 edge weights.

Functions
---------
- (ti,tj,tv,nverts) = mst_prim(A::MatrixNetwork,full::Bool,u::Int64)
"""
SCOMPONENTS
-----------
compute the strongly connected components of a graph

ci=scomponents(A) returns an index for the component number of every 
vertex in the graph A.  The total number of components is maximum(ci).
If the input is undirected, then this algorithm outputs just the 
connected components.  Otherwise, it output the strongly connected components.

The implementation is from Tarjan's 1972 paper: Depth-first search and 
linear graph algorithms. In SIAM's Journal of Computing, 1972, 1, 
pp.146-160.

Functions
---------
- cc = scomponents(A::MatrixNetwork)
- cc = scomponents{T}(A::SparseMatrixCSC{T,Int64}) 
- sci = strong_components_map(A::MatrixNetwork)
- sci = strong_components_map{T}(A::SparseMatrixCSC{T,Int64})
- sc_rich = enrich(cc::Strong_components_output) # check ?enrich for more

Example
-------
A = load_matrix_network("cores_example")
cc = scomponents(A)
scomponents(A).number
scomponents(A).sizes
scomponents(A).map
strong_components_map(A)     # if you just want the map
enrich(scomponents(A)) # produce additional enriched output

Can work on ei,ej\n
ei = [1;2;3]
ej = [2;4;1]
scomponents(ei,ej)

Can work on sparse matrix A\n
A = sprand(5,5,0.5)
scomponents(A)
"""

:scomponents

###########################
##    Type Definitions    #
###########################
type Strong_components_output
    map::Vector{Int64}    # The indicator map
    sizes::Vector{Int64}  # Array of sizes corresponding to map
    number::Int64         # Int64
    A::MatrixNetwork      # MatrixNetwork
end

type Strong_components_rich_output
    reduction_matrix::SparseMatrixCSC{Int64,Int64} # the reduction matrix (restriction matrix)
    transitive_order::SparseMatrixCSC{Int64,Int64}
    transitive_map::SparseMatrixCSC{Int64,Int64}
end

####################
##    Functions    #
####################


"""
Return information on the strongly connected components of a graph that 
is the minimum required computation.
Example:
``strong_components_map(MatrixNetwork(sprand(5,4,0.1)))``	
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

########################################################
##    Conversion Functions for strong_components_map   #
########################################################
# CSC:
strong_components_map{T}(A::SparseMatrixCSC{T,Int64}) = 
                                        strong_components_map(MatrixNetwork(A))
# Triplet:
strong_components_map(ei::Vector{Int64},ej::Vector{Int64}) = 
                                        strong_components_map(MatrixNetwork(ei,ej))
# CSR
strong_components_map{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64) = 
                                        strong_components_map(MatrixNetwork(n,rp,ci,vals))

######################
##    scomponents    #
######################

function scomponents(A::MatrixNetwork)
    map = strong_components_map(A)
	number = maximum(map)
	sizes = zeros(Int64,number);
    for i = 1:length(map)
        sizes[map[i]] += 1;
    end
    
    return Strong_components_output(map, sizes, number, A)
end

###############################
##    Conversion Functions    #
###############################

# CSC:
scomponents{T}(A::SparseMatrixCSC{T,Int64}) = scomponents(MatrixNetwork(A))
# Triplet:
scomponents(ei::Vector{Int64},ej::Vector{Int64}) = scomponents(MatrixNetwork(ei,ej))
# CSR
scomponents{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64) = scomponents(MatrixNetwork(n,rp,ci,vals))

""" 
This function adds the following helpers variables
* reduction_matrix - a reduction matrix to project down to the component graph
* transitive_order - a transitive ordering of the components
* transitive_map - a map to components that respects the transitive ordering
* largest
"""
function enrich(rval::Strong_components_output)
    ci = rval.map
    sizes = rval.sizes
    ncomp = maximum(ci)
    R = sparse(collect(1:rval.A.n),ci,1,rval.A.n,ncomp)
    A = SparseMatrixCSC(rval.A.n,rval.A.n,rval.A.rp,rval.A.ci,rval.A.vals)
    CG = R'*A'*R
    return Strong_components_rich_output(R,A,CG)
end





