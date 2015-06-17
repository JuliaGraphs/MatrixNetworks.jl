
# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation
# TODO: add more examples
# TODO support struct


"""
Example
-------

using MAT

file_path = Pkg.dir("MatrixNetworks/data/dfs_example.mat")

file = matopen(file_path)

A = read(file,"A")

close(file)

dfs(MatrixNetwork(A),1)

DFS compute depth first search distances and returns the distance (d), the discover (dt),
the finish time(ft), and the predecessor array (pred) in the tuple (d,dt,ft, pred).

pred[i] = 0 if vertex i is in a component not reachable from u and i != u.




"""

function dfs(A::MatrixNetwork,u::Int64)
    return dfs(A,u,0,0);
end

function dfs(A::MatrixNetwork,u::Int64,full::Int64,target::Int64)

    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    d=-1*ones(Int64,n)
    dt=-1*ones(Int64,n)
    ft=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    rs=zeros(Int64,2*n)
    rss=0 # recursion stack holds two nums (v,ri)
    
    #start dfs at u
    t=0
    targethit=0
    for i=1:n
        if i==1
            v=u
        else
            v=mod(u+i-1,n)+1
            if d[v]>0
                continue
            end
        end
        d[v]=0
        dt[v]=t
        t=t+1
        ri=rp[v]
        rss=rss+1
        rs[2*rss-1]=v
        rs[2*rss]=ri # add v to the stack
        
        while rss>0
            v=rs[2*rss-1]
            ri=rs[2*rss]
            rss=rss-1 # pop v from the stack
            if v==target || targethit == 1
                ri=rp[v+1]
                targethit=1 # end the algorithm if v is the target
            end
            while ri<rp[v+1]
                w=ci[ri]
                ri=ri+1
                if d[w]<0
                    d[w]=d[v]+1
                    pred[w]=v
                    rss=rss+1
                    rs[2*rss-1]=v
                    rs[2*rss]=ri # add v to the stack
                    v=w
                    ri=rp[w]
                    dt[v]=t
                    t=t+1
                    continue #discover a new vertex!
                end
            end
            ft[v]=t
            t=t+1 # finish with v
        end
        if full == 0
            break
        end
    end
    return (d,dt,ft,pred)
end
    
