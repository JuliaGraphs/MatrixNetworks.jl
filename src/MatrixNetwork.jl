# create type MatrixNetwork
type MatrixNetwork{T}
    n::Int64 # number of columns/rows
    rp::Vector{Int64} # row pointers
    ci::Vector{Int64} # column indices
    vals::Vector{T} # corresponding values
end

import Base.sparse

"""
Return back an adjacency matrix representation
of the transpose. This
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
