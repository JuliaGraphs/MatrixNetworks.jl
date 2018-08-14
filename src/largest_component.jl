"""
LARGEST_COMPONENT
-----------------
    return the largest connected component of A.
    Acc = largest_component(A) returns the largest connected component
    of the graph A.  If A is directed, this returns the largest
    strongly connected component.
    
    Acc = largest_component(A,true) returns the largest connected piece of
    a directed graph where connectivity is undirected.  Algorithmically,
    this takes A, drops the directions, then components the largest component
    and returns just this piece of the original _directed_ network.  So the
    output Acc is directed in this case.
    
    (Acc,p) = largest_component(A) also returns a logical vector
    indicating which vertices in A were chosen.

Functions
---------
- (Acc,p) = largest_component{T}(A::SparseMatrixCSC{T,Int64})
- (Acc,p) = largest_component{T}(A::SparseMatrixCSC{T,Int64},sym::Bool)

Example
-------
~~~
A = load_matrix_network("dfs_example")
(Acc,p) = largest_component(A)
~~~
"""
function largest_component(A::SparseMatrixCSC{T,Int64}) where T
    return largest_component(A,false)
end

function largest_component(A::SparseMatrixCSC{T,Int64},sym::Bool) where T
    if sym
        # As = A|A' until Julia implements this on sparse matrices, use this:
        As = max.(A,A')
        cc = scomponents(As)
    else
        cc = scomponents(A)
    end
    cind = argmax(cc.sizes)
    p = cc.map .== cind
    # no logical indexing so:
    idx = findall(p)
    Acc = A[idx,idx]
    return (Acc,p)
end

    
