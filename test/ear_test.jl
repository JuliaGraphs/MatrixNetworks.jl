using MatrixNetworks

"""

Ear Decomposition
-----------------

An ear decomposition of an undirected graph G is a partition of its set of edges into a sequence of ears, such that the one or two endpoints of each ear 
belong to earlier ears in the sequence and such that the internal vertices of each ear do not belong to any earlier ear.

An open ear decomposition or a proper ear decomposition is an ear decomposition in which the two endpoints of each ear after the first are distinct from each other.

The code is based on the algorithm given by Schmidt et.al in  "A simple test on 2-vertex- and 2-edge-connectivity" (Information Processing Letters Journal'
2013).

This approach is an important subroutine in many graph problems like Planarity Testing, Identifying the Factor critical graphs, recognizing series parallel 
graphs, identifying Tri connected components etc.

Input: An undirected graph with i--->j and j--->i present (NO self loops!)

Observed:

      1) Slower on random graphs generated using GT-graph Suite ( Erdos Reyni model - with 20K vertices and a sparse quotient of 10 ) by a factor 5-8%

TODO: 1) Add a function to test for bi-connectivity and 2-edge connectivity.
      2) Compare the Performance with C ( SNAP DATA set )

"""


I = [1,2,2,3,3,4,1,4,3,5,5,6,6,7,7,8,6,8,5,9,9,10,5,10,1,3,2,5]
J = [2,1,3,2,4,3,4,1,5,3,6,5,7,6,8,7,8,6,9,5,10,9,10,5,3,1,5,2]
V = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

A = sparse(I,J,V,10,10)

function ear(A::MatrixNetwork,u::Int64,full::Int64,target::Int64)

    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    d=-1*ones(Int64,n)
    dt=-1*ones(Int64,n)
    rev_dt= -1*ones(Int64,n)
    ft=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    rs=zeros(Int64,2*n)
    rss=0 # recursion stack holds two nums (v,ri)
    be = Any[] #back edges array   
    ears = zeros(Int64,2 * length(A.ci))
    ear_number = zeros(Int64,length(A.ci))
    ear_index = 0
    ear_visited = zeros(Int64, A.n)

    for i=1:n
        push!(be, [])
    end

    #start dfs at u
    t=1
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
        d[v]=1
        dt[v]=t
        rev_dt[t]=v;
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
		    rev_dt[t]=v
                    t=t+1
                #    continue #discover a new vertex!
                elseif (pred[v]!=w)
		    if(dt[w]< dt[v])
		       push!(be[w], v)
		    end
		end	
            end
        end
        if full == 0
            break
        end
    end
    cn=1
    for i=1:n
	startVertex = rev_dt[i]
	ear_visited[startVertex] = 1
	while !isempty(be[startVertex])  # A disjoint Path (or edge) is identified.
	    nextVertex = pop!(be[startVertex])
            ear_index += 1
            ears[2*ear_index-1] = startVertex
            ears[2*ear_index] = nextVertex
            ear_number[ear_index] = cn
            while (ear_visited[nextVertex] != 1) #Identify the ear (Chain) and number it
                ear_visited[nextVertex] = 1
                startVertex_n = nextVertex
                nextVertex = pred[startVertex_n]    
                ear_index += 1
                ears[2*ear_index-1] = startVertex_n # Store the ear (Chain). 
                ears[2*ear_index] = nextVertex
                ear_number[ear_index] = cn  # Provide numbering to that Edge (#ear)
            end 
            cn += 1	
	end
    end
    print(ears)
    return (ears,ear_number)
end


function dfs_test()
#print(A);
B = MatrixNetwork(A);
#print(B) 
ear(B, 1,1,0)
end

dfs_test();
