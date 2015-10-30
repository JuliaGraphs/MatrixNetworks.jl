"""
CSRTOSPARSE
-----------
    convert a matrix from compressed sparse row to a sparse matrix A.
    It returns the arrays that feed the sparse function in julia.

Functions
---------
- (nzi,nzj,nzv) = csr_to_sparse{T}(rp::Vector{Int64},ci::Vector{Int64},ai::Vector{T})
- (nzi,nzj,nzv) = csr_to_sparse{T}(rp::Vector{Int64},ci::Vector{Int64},ai::Vector{T},nrows::Int64)

Example
-------
~~~
i = [1;2;3]
j = [3;4;4]
v = [8;9;10]
(rp,ci,ai,m) = sparse_to_csr(i,j,v)
(nzi,nzj,nzv) = csr_to_sparse(rp,ci,ai)
A = sparse(nzi,nzj,nzv,length(rp)-1,maximum(ci))
B = csr_to_sparse_matrix(rp,ci,ai)
isequal(A,B)
~~~
"""

##################
#	Functions    #
##################
function csr_to_sparse{T}(rp::Vector{Int64},ci::Vector{Int64},ai::Vector{T})
    nrows = length(rp)-1
    ncols = length(ci)
    nzi = zeros(Int64,ncols)
    for i=1:nrows
        for j=rp[i]:rp[i+1]-1
            nzi[j] = i
        end
    end
    nzj = ci
    nzv = ai
    return (nzi,nzj,nzv)
end 

function csr_to_sparse_matrix{T}(rp::Vector{Int64},ci::Vector{Int64},ai::Vector{T})
    (i,j,k) = csr_to_sparse(rp,ci,ai)
    A = sparse(i,j,k)
    return A
end

function csr_to_sparse_matrix{T}(rp::Vector{Int64},ci::Vector{Int64},
                                            ai::Vector{T},nrows::Int64,ncols::Int64)
    (i,j,k) = csr_to_sparse(rp,ci,ai);
    A = sparse(i,j,k,nrows,ncols)
    return A
end

