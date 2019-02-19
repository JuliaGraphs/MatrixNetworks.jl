#TODO: more detailed documentation?
"""
BIPARTITE MATCHING
------------------
    return a maximum weight bipartite matching of a graph

Functions
---------
- Matching_Setup = bipartite_matching_setup{T}(A::SparseMatrixCSC{T,Int64})
- Matching_Setup = bipartite_matching_setup{T}(x::Vector{T},ei::Vector{Int64},ej::Vector{Int64},m::Int64,n::Int64)
- Matching_Output = bipartite_matching_primal_dual{T}(rp::Vector{Int64}, ci::Vector{Int64},ai::Vector{T}, m::Int64, n::Int64)
- Matching_Output = bipartite_matching{T}(A::SparseMatrixCSC{T,Int64})
- Matching_Output = bipartite_matching{T}(w::Vector{T},ei::Vector{Int64},ej::Vector{Int64},m::Int64,n::Int64)
- Matching_Output = bipartite_matching{T}(w::Vector{T},ei::Vector{Int64},ej::Vector{Int64})
- ind = bipartite_matching_indicator{T}(w::Vector{T},ei::Vector{Int64},ej::Vector{Int64})
- (m1,m2) = edge_list(M_output::Matching_output)
- S = create_sparse(M_output::Matching_output)\n
You can check the documentation of each of the output modifiers functions separately.

Example
-------
~~~
W = sprand(10,8,0.5)
bipartite_matching(W)
ei = [1;2;3]
ej = [3;2;4]
Matching_Output = bipartite_matching([10;12;13],ei,ej)
Matching_Output.weight
Matching_Output.cardinality
Matching_Output.match
S = create_sparse(bipartite_matching(W)) # get the sparse matrix
(m1,m2) = edge_list(bipartite_matching(W)) # get the edgelist
~~~
"""
function bipartite_matching end

mutable struct Matching_setup
    rp::Array{Int64,1}
    ci::Array{Int64,1}
    ai::Array{Float64,1}
    tripi::Array{Int64,1}
    m::Int64
    n::Int64
end

mutable struct Matching_output
    m::Int64
    n::Int64
    weight::Float64
    cardinality::Int64
    match::Array{Int64,1}
end

######################
#   setup  funtions  #
######################

function bipartite_matching_setup_phase1(A::SparseMatrixCSC{T,Int64}) where T
    (nzi,nzj,nzv) = findnz(A)
    return (nzi,nzj,nzv)
end

function bipartite_matching_setup(A::SparseMatrixCSC{T,Int64}) where T
    (nzi,nzj,nzv) = bipartite_matching_setup_phase1(A)
    nedges = length(nzi)
    (m,n)=size(A)
    
    rp = ones(Int64,m+1) # csr matrix with extra edges
    ci = zeros(Int64,nedges+m)
    ai = zeros(Float64,nedges+m)
    
    
    rp[1]=0
    for i=1:nedges
        rp[nzi[i]+1]=rp[nzi[i]+1]+1
    end
    rp=cumsum(rp) 
    
    for i=1:nedges
        ai[rp[nzi[i]]+1]=nzv[i]
        ci[rp[nzi[i]]+1]=nzj[i]
        rp[nzi[i]]=rp[nzi[i]]+1
    end
    
    for i=1:m # add the extra edges
        ai[rp[i]+1]=0
        ci[rp[i]+1]=n+i
        rp[i]=rp[i]+1
    end
    
    # restore the row pointer array
    for i=m:-1:1
        rp[i+1]=rp[i]
    end
    rp[1]=0
    rp=rp.+1
    
    #check for duplicates in the data
    colind = zeros(Int64,m+n)
    for i=1:m
        for rpi=rp[i]:rp[i+1]-1
            if colind[ci[rpi]] == 1
                error("bipartite_matching:duplicateEdge")
            end
        colind[ci[rpi]]=1
        end
    
        for rpi=rp[i]:rp[i+1]-1
            colind[ci[rpi]]=0
        end # reset indicator
    end
    
    M_setup = Matching_setup(rp,ci,ai,[],m,n)
    return M_setup
end


function bipartite_matching_setup(x::Vector{T},ei::Vector{Int64},
                                     ej::Vector{Int64},m::Int64,n::Int64) where T
    (nzi,nzj,nzv) = (ei,ej,x)
    nedges = length(nzi)
    
    rp = ones(Int64,m+1) # csr matrix with extra edges
    ci = zeros(Int64,nedges+m)
    ai = zeros(Float64,nedges+m)
    tripi = zeros(Int64,nedges+m)
    # 1. build csr representation with a set of extra edges from vertex i to
    # vertex m+i
    
    rp[1]=0
    for i=1:nedges
        rp[nzi[i]+1]=rp[nzi[i]+1]+1
    end
    
    rp=cumsum(rp)
    for i=1:nedges
        tripi[rp[nzi[i]]+1]=i
        ai[rp[nzi[i]]+1]=nzv[i]
        ci[rp[nzi[i]]+1]=nzj[i]
        rp[nzi[i]]=rp[nzi[i]]+1
    end
    for i=1:m # add the extra edges
        tripi[rp[i]+1]=-1
        ai[rp[i]+1]=0
        ci[rp[i]+1]=n+i
        rp[i]=rp[i]+1
    end
    
    # restore the row pointer array
    for i=m:-1:1
        rp[i+1]=rp[i]
    end
    
    rp[1]=0
    rp=rp.+1
    
    #check for duplicates in the data
    colind = zeros(Int64,m+n)
    
    for i=1:m
        for rpi=rp[i]:rp[i+1]-1
            if colind[ci[rpi]] == 1
                error("bipartite_matching:duplicateEdge")
            end
        colind[ci[rpi]]=1
        end
    
        for rpi=rp[i]:rp[i+1]-1
            colind[ci[rpi]]=0 
        end # reset indicator
    end
    M_setup = Matching_setup(rp,ci,ai,tripi,m,n)
    return M_setup
end

##################
#   primal-dual  #
##################

function bipartite_matching_primal_dual(rp::Vector{Int64}, ci::Vector{Int64}, 
                    ai::Vector{T}, m::Int64, n::Int64) where T
    
    # variables used for the primal-dual algorithm
    # normalize ai values # updated on 2-19-2019
    ai ./= maximum(abs.(ai))
    alpha=zeros(Float64,m)
    bt=zeros(Float64,m+n)#beta
    queue=zeros(Int64,m)
    t=zeros(Int64,m+n)
    match1=zeros(Int64,m)
    match2=zeros(Int64,m+n)
    tmod = zeros(Int64,m+n)
    ntmod=0
    
    # initialize the primal and dual variables
    for i=1:m
        for rpi=rp[i]:rp[i+1]-1
            if ai[rpi] > alpha[i]
               alpha[i]=ai[rpi]
            end
        end
    end
    
    # dual variables (bt) are initialized to 0 already
    # match1 and match2 are both 0, which indicates no matches
    
    i=1
    while i<=m
        for j=1:ntmod
            t[tmod[j]]=0
        end
        ntmod=0
        # add i to the stack
        head=1
        tail=1
        queue[head]=i
        while head <= tail && match1[i]==0
            k=queue[head]
            for rpi=rp[k]:rp[k+1]-1
                j = ci[rpi]
                if ai[rpi] < alpha[k] + bt[j] - 1e-8
                    continue
                end # skip if tight
                if t[j]==0
                    tail=tail+1
                    if tail <= m
                        queue[tail]=match2[j]
                    end
                    t[j]=k
                    ntmod=ntmod+1
                    tmod[ntmod]=j
                    if match2[j]<1
                        while j>0
                            match2[j]=t[j]
                            k=t[j]
                            temp=match1[k]
                            match1[k]=j
                            j=temp
                        end
                        break
                    end
                end
            end
            head=head+1
        end
        if match1[i] < 1
            theta=Inf
            for j=1:head-1
                t1=queue[j]
                for rpi=rp[t1]:rp[t1+1]-1
                    t2=ci[rpi]
                    if t[t2] == 0 && alpha[t1] + bt[t2] - ai[rpi] < theta
                        theta = alpha[t1] + bt[t2] - ai[rpi]
                    end
                end
            end
            for j=1:head-1
                alpha[queue[j]] -= theta
            end
            for j=1:ntmod
                bt[tmod[j]] += theta
            end
            continue
        end
        i=i+1
    end
    val=0
    for i=1:m
        for rpi=rp[i]:rp[i+1]-1
            if ci[rpi]==match1[i]
                val=val+ai[rpi]
            end
        end
    end
    noute = 0
    for i=1:m
        if match1[i]<=n
            noute=noute+1
        end
    end

    M_output = Matching_output(m,n,val,noute,match1)
    return M_output
end

function bipartite_matching_primal_dual(M_setup::Matching_setup)
    return bipartite_matching_primal_dual(M_setup.rp, M_setup.ci, M_setup.ai,
                                          M_setup.m, M_setup.n)
end

####################
##    Functions    #
####################

function bipartite_matching(A::SparseMatrixCSC{T,Int64}) where T
    return bipartite_matching_primal_dual(bipartite_matching_setup(A))
end


function bipartite_matching(w::Vector{T},ei::Vector{Int64},
                                     ej::Vector{Int64},m::Int64,n::Int64) where T
    return bipartite_matching_primal_dual(bipartite_matching_setup(w,ei,ej,m,n))
end

function bipartite_matching(w::Vector{T},ei::Vector{Int64},
                                     ej::Vector{Int64}) where T
    return bipartite_matching_primal_dual(bipartite_matching_setup(
        w,ei,ej,maximum(ei),maximum(ej)))
end

####################
##    Indicator    #
####################

"""
Returns the matching indicator of a matrix stored in triplet format
Example:
bipartite_matching_indicator([10;12;13],[1;2;3],[3;2;4])
"""
function bipartite_matching_indicator(w::Vector{T},ei::Vector{Int64},
                                     ej::Vector{Int64}) where T
    M_setup = bipartite_matching_setup(w, ei, ei, maximum(ei), maximum(ej))
    M_out = bipartite_matching_primal_dual(M_setup.rp, M_setup.ci, M_setup.ai,
                                           M_setup.m, M_setup.n)
    return edge_indicator(M_out,ei, ej)
end


#########################
#   Output Modifiers    #
#########################
"""
Returns the edge list of a matching output
Example:
M_out = bipartite_matching([10;12;13],[1;2;3],[3;2;4])
edge_list(M_out)
"""
function edge_list(M_output::Matching_output)
    m1=zeros(Int64,M_output.cardinality)
    m2=zeros(Int64,M_output.cardinality)
    noute=1
    for i=1:M_output.m
        if M_output.match[i]<=M_output.n
            m1[noute]=i
            m2[noute]=M_output.match[i]
            noute=noute+1
        end
    end
    return (m1,m2)
end

##########
"""
Creates and returns a sparse matrix that represents the outputed matching
Example:
M_out = bipartite_matching([10;12;13],[1;2;3],[3;2;4])
create_sparse(M_out)
"""
function create_sparse(M_output::Matching_output)
    (ei,ej) = edge_list(M_output)
    A = sparse(ei,ej,1,M_output.m,M_output.n)
    return A
end

##########
function edge_indicator(M_output::Matching_output, ei::Vector, ej::Vector)
    @assert length(ei) == length(ej)
    ind = zeros(Int64,length(ei))
    for i=1:length(ei)
        if M_output.match[ei[i]] == ej[i]
            ind[i] = 1
        end
    end
    return ind
end
