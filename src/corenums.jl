"""
CORENUMS
---------
    compute the core number for each vertex in the graph and returns the core
    numbers for each vertex of the graph A along with the removal order of the vertex in the 
    tuple (d,rt). This function works on directed graphs but gives the in-degree core number.
    To get the out-degree core numbers call corenums(A')

Functions
---------
- (d,rt) = corenums(A::MatrixNetwork)
- (d,rt) = corenums{T}(A::SparseMatrixCSC{T,Int64})
- (d,rt) = corenums(ei::Vector{Int64},ej::Vector{Int64})

Example
-------
~~~
A = load_matrix_network("cores_example")
(d,rt) = corenums(A)
~~~
"""
function corenums(A::MatrixNetwork)
    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    nz=length(ci)

    d=zeros(Int64,n)
    maxd=0
    rt=zeros(Int64,n)
    for k=1:nz
        newd=d[ci[k]]+1
        d[ci[k]]=newd
        if newd>maxd
            maxd=newd
        end
    end
    
    # compute the bucket sort
    dp=zeros(Int64,maxd+2)
    vs=zeros(Int64,n)
    vi=zeros(Int64,n)# degree position, vertices
    for i=1:n
        dp[d[i]+2]=dp[d[i]+2]+1
    end # plus 2 because degrees start at 0
    
    dp=cumsum(dp)
    dp=dp.+1
    for i=1:n
        vs[dp[d[i]+1]]=i
        vi[i]=dp[d[i]+1]
        dp[d[i]+1]=dp[d[i]+1]+1
    end
    
    for i=maxd:-1:1
        dp[i+1]=dp[i]
    end
    
    # start the algorithm
    t=1
    for i=1:n
        v = vs[i]
        dv = d[v]
        rt[v]=t
        t=t+1
        
        for rpi=rp[v]:rp[v+1]-1
            w=ci[rpi]
            dw=d[w]
            
            if dw<=dv # we already removed w
            else #need to remove edge (v,w), which decreases d(w)
                # swap w with the vertex at the head of its degree
                pw=vi[w] # get the position of w
                px=dp[dw+1] #get the pos of the vertex at the head of dw list
                x=vs[px]
                #swap w, x
                vs[pw]=x
                vs[px]=w
                vi[w]=px
                vi[x]=pw
                
                # decrement the degree of w and increment the start of dw
                d[w]=dw-1
                dp[dw+1]=px+1
            end
        end
    end
    return (d,rt)
end

############################
### Additional functions ###
############################

## CSC:
function corenums(A::SparseMatrixCSC{T,Int64}) where T
    return corenums(MatrixNetwork(A))
end

## Triplet Format:
function corenums(ei::Vector{Int64},ej::Vector{Int64})
    return corenums(MatrixNetwork(ei,ej))
end

## CSR sparse matrices:
function corenums(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64) where T
    return corenums(MatrixNetwork(n,rp,ci,vals))
end