# Use the new docile convention
# http://docilejl.readthedocs.org/en/latest/syntax/
# TODO: more testing and check documentation
# TODO: add more examples
# TODO support struct input
# TODO support other similarity functions
# TODO allow option for diagonal similarity too
# TODO allow option for triplet output


"""
Example
-------

S = corenums(A)

"""


"""
cosineknn compute the k-nearest neighbors similarity metric between the
vertices of A or the upper half of a bipartite graph A
"""
:cosineknn

function cosineknn(A::SparseMatrixCSC{Float64,Int64},k::Int64)
end



# (rp,ci,ai) = sparse_to_csr(A)
# check = 1
# (rpt,cit,ait) = (A.colptr,A.rowval,A.nzval)
# (n,m) = size(A)

