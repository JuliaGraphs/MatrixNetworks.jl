# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation

"""
Example
-------

i = [1;2;3];
j = [3;4;4];
v = [8;9;10];
(rp,ci,ai,m) = sparse_to_csr(i,j,v)
(nzi,nzj,nzv) = csr_to_sparse(rp,ci,ai)

# To construct the sparse matrix in julia:
A = sparse(nzi,nzj,nzv,length(rp)-1,maximum(ci))

"""


##################
#	Functions    #
##################

"""
csr_to_sparse converts a matrix from compressed sparse row to a sparse matrix A
it returns the arrays that feed the sparse function in julia.
"""
:csr_to_sparse
function csr_to_sparse{T}(rp::Array{Int64,1},ci::Array{Int64,1},ai::Array{T,1},varargin...)
    if length(varargin)==0
        nrows = length(rp)-1;
    else
        nrows = varargin[1];
    end

    ncols = length(ci);
    nzi = zeros(Int64,ncols);
    for i=1:nrows
        for j=rp[i]:rp[i+1]-1
            nzi[j] = i;
        end
    end
    nzj = ci;
    nzv = ai;
    return (nzi,nzj,nzv)
end 