"""
SPARSETOCSR
-----------
    convert sparse matrix into compressed row storage arrays
    and returns the row pointer (rp), column index (ci) and value index (ai) arrays 
    of a compressed sparse representation of the matrix A (or its triplet format storage)

Functions
---------
- (rp,ci,ai,m) = sparse_to_csr{T}(A::SparseMatrixCSC{T,Int64})
- (rp,ci,ai,m) = sparse_to_csr{T}(nzi::Array{Int64,1},nzj::Array{Int64,1},nzv::Array{T,1}) 

Example
-------
~~~
i = [1;2;3]
j = [3;4;4]
v = [8;9;10]
(rp,ci,ai,m) = sparse_to_csr(i,j,v)
~~~
"""
function sparse_to_csr(A::SparseMatrixCSC{T,Int64}) where T
    At = copy(A');
    return (At.colptr,At.rowval,At.nzval,At.m);
end

function sparse_to_csr(nzi::Array{Int64,1},nzj::Array{Int64,1},
                         nzv::Array{T,1}) where T 
    At = sparse(nzj,nzi,nzv)
    return (At.colptr,At.rowval,At.nzval,At.m)
end


import SparseArrays: findnz   #import function to extend to MatrixNetwork type
"""
findnz
======
Returns the uncompressed non-zeros of the CSR representation of a MatrixNetwork.
Js and Vs are copied from the CSR format, and Is is computed. 

Output
------
- Is::Array{Integer,1}: the row indices.
- Js::Array{Integer,1}: the columb indices.
- Vs::Array{Integer,1}: the non-zero values.

Example
-------
Is, Js, Vs = findnz(A)
"""
function findnz(A::MatrixNetwork{T}) where T

    Is = zeros(Int64,length(A.vals))
    Js = copy(A.ci)
    Vs = copy(A.vals)

    idx = 1
    for i in 1:(length(A.rp)-1)
        row_nz = A.rp[i+1] - A.rp[i]
        for _ in 1:row_nz
            Is[idx] = i
            idx += 1 
        end
    end

    return Is,Js,Vs
end