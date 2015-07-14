"""
Example: 
A = load_matrix_network("dfs_example")
(Acc,p) = largest_component(A)
"""

function largest_component{T}(A::SparseMatrixCSC{T,Int64})
    return largest_component(A,false)
end

function largest_component{T}(A::SparseMatrixCSC{T,Int64},sym::Bool)
    if sym
        As = A|A'
        cc = scomponents(As)
    else
        cc = scomponents(A)
    end
    cind = indmax(cc.sizes)
    p = cc.map .== cind
    # no logical indexing so:
    idx = find(p)
    Acc = A(idx,idx)
    return (Acc,p)
end

    
