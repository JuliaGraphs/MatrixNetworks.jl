
# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation
# TODO: add more examples
# TODO support struct

"""
Example
-------

using MAT

file_path = Pkg.dir("MatrixNetworks/data/cores_example.mat")
  
file = matopen(file_path)

A = read(file,"A")

close(file)

(d,rt) = corenums(MatrixNetwork(A))

corenums compute the core number for each vertex in the graph and returns the core
numbers for each vertex of the graph A along with the removal order of the vertex in the 
tuple (d,rt)
This function works on directed graphs but gives the in-degree core number.
To get the out-degree core numbers call corenums(A')
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
    dp=dp+1
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
