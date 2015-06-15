# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation

"""
Example
-------

i = [1;2;3];
j = [3;4;4];
v = [8;9;10];
(rp,ci,ai,m) = sparse_to_csr(i,j,v,4)

"""



#################
#   Functions   #
#################
"""
sparse_to_csr converts a sparse matrix into compressed row storage arrays
and returns the row pointer (rp), column index (ci) and value index (ai) arrays 
of a compressed sparse representation of the matrix A (or its triplet format storage)
"""
:sparse_to_csr
function sparse_to_csr{T}(A::SparseMatrixCSC{T,Int64})
    At = A';
    return (At.colptr,At.rowval,At.nzval,At.m);
end

function sparse_to_csr{T}(nzi::Array{Int64,1},nzj::Array{Int64,1},
                            nzv::Array{T,1},varargin...) 
    At = sparse(nzj,nzi,nzv)
    return (At.colptr,At.rowval,At.nzval,At.m)
end
