using Compat 

# create type MatrixNetwork
type MatrixNetwork{T}
    n::Int64 # number of columns/rows
    rp::Vector{Int64} # row pointers
    ci::Vector{Int64} # column indices
    vals::Vector{T} # corresponding values
end

function MatrixNetwork{T}(A::SparseMatrixCSC{T,Int64})
    At = A'
    return MatrixNetwork(size(At,2),At.colptr,At.rowval,At.nzval)
end

MatrixNetwork(ei::Vector{Int64},ej::Vector{Int64}) = 
    MatrixNetwork(ei,ej,max(maximum(ei),maximum(ej)))

function MatrixNetwork(ei::Vector{Int64},ej::Vector{Int64},n::Int64)
    At = sparse(ej,ei,true,n,n);
    return MatrixNetwork(size(At,2),At.colptr,At.rowval,At.nzval)
end


function _matrix_network_direct{T}(A::SparseMatrixCSC{T,Int64})
    return MatrixNetwork(size(A,2),A.colptr,A.rowval,A.nzval)
end

function _matrix_network_direct{T}(A::SparseMatrixCSC{T,Int64},v)
    nzval = ones(typeof(v),length(A.nzval))
    return MatrixNetwork(size(A,2),A.colptr,A.rowval,nzval)
end

import Base.sparse, Base.size

"""
Return back an adjacency matrix representation
of the transpose. This requires no work. 
"""
function sparse_transpose{T}(A::MatrixNetwork{T})
    return SparseMatrixCSC(A.n,A.n,A.rp,A.ci,A.vals)
end

"""
Return back an adjacency matrix representation
of the current MatrixNetwork
"""
function sparse{T}(A::MatrixNetwork{T})
    return sparse_transpose(A)'
end

function size(A::MatrixNetwork)
    return (A.n,A.n)
end

function size(A::MatrixNetwork, dim::Integer)
    if dim == 1 || dim == 2
        return A.n
    elseif dim > 2
        return 1
    else
        throw(DomainError())
    end
end

"""
`is_undirected`
===============

Check the matrix associated with a matrix network
for symmetry. 

Input
-----
- `A`: a matrix network

Returns
-------
- `bool` with true indicating the network is undirected
    and the matrix is symmetric
"""    
function is_undirected end
    
function is_undirected(A::MatrixNetwork)
   M = sparse_transpose(A)
   return issymmetric(M) 
end

function is_undirected(A::SparseMatrixCSC)
   return issymmetric(A) 
end


"""
`is_connected`
==============

Check the matrix associated with a matrix network
for (strong) connectivity  

Usage
-----
- `is_connected(A)`

Input
-----
- `A`: a `MatrixNetwork` or `SparseMatrixCSC` class

Returns
-------
- `bool` with true indicating the matrix is strongly connected
and false indicating 
"""    
function is_connected end

function is_connected(A::MatrixNetwork)
    return maximum(scomponents(A).map) == 1
end

function is_connected(A::SparseMatrixCSC)
    return maximum(scomponents(A).map) == 1
end