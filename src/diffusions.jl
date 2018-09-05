
## Todo

import LinearAlgebra.checksquare
import Base.eltype
import Base.length
import Base.ndims
import Base.*
import LinearAlgebra.mul!
import Base.size
import KahanSummation

"""
A simple function that doesn't report any output
"""
function _noiterfunc(iter::Int, x::Vector{T}) where T
end

"""
`_applyv!`
--------

The fully generic function to compute
`x <- alpha*x + gamma*v`

Functions
---------
- `_applyv!(x, v, alpha, gamma)` 

Example
-------
Internal function
"""
function _applyv! end

function _applyv!(x::Vector{T}, v::SparseVector{T},
                 alpha::T, gamma::T) where T
    @simd for i in 1:length(x)
        @inbounds x[i] *= alpha
    end
    vvals = v.nzval
    vrows = v.nzind
    @simd for j in 1:length(vrows)
        @inbounds x[vrows[j]] += gamma*vvals[j]
    end
end
 # version 0.5.0-dev

function _applyv!(x::Vector{T}, v, alpha::T, gamma::T) where T
    x *= alpha
    x += gamma*v
end

function _applyv!(x::Vector{T}, v::T, alpha::T, gamma::T) where T
    gv = gamma*v
    @simd for i in 1:length(x)
        @inbounds x[i] = alpha*x[i] + gv
    end
end

function _applyv!(x::Vector{T}, v::Vector{T}, alpha::T, gamma::T) where T
    @simd for i in 1:length(x)
        @inbounds x[i] = alpha*x[i] + gamma*v[i]
    end
end

function _applyv!(x::Vector{T}, v::SparseMatrixCSC{T,Int},
                 alpha::T, gamma::T) where T
    @simd for i in 1:length(x)
        @inbounds x[i] *= alpha
    end
    vvals = nonzeros(v)
    vrows = rowvals(v)
    @simd for j in nzrange(v,1)
        @inbounds x[vrows[j]] += gamma*vvals[j]
    end
end

"""
Convert a sparse representation into a dense vector.
"""
function _densevec(I::Vector{Int}, V::Vector{T}, n::Int) where T
    v = zeros(T,n)
    for ind in 1:length(I)
        v[ind] = V[ind]
    end
    return v
end
function _densevec(d::Dict{Int,T}, n::Int) where T
    v = zeros(T,n)
    for (ind,val) in d
        v[ind] = val
    end
    return v
end

"""
Return true if any element is negative.
"""
function _check_negative(v)
    for val in v
        if val < 0
            return true
        end        
    end 
    return false
end



"""
`pagerank_power!`
----------------

**Expert interface**

This function computes the strongly personalized PageRank
vector of a column sub-stochastic matrix P.

This is a powerful internal function and you most likely don't want 
to use it unless you know what you are doing.  

This call is duck-typed so that it can work with arbitrary input
matrix types as well as arbitrary vectors v. This enables it to
be efficient in the case of sparse or dense vectors v.

The solution returned will satisfy the strongly-personalized
PageRank equation to the given tolerance accuracy in the 1-norm.
Formally, it'll 

This call allocates no extra memory and only uses the memory
included with the vectors x and y. 

Functions
---------
- `pagerank_power!{T<Float}!(x::Vector{T}, y::Vector{T}, 
        P, alpha::T, v, tol::T, maxiter::Int, iterfunc::Function)` 
                
Inputs
------
- `x`: a vector (of length n, where there are n nodes of the graph)
  that will store 
- `y`: a vector that is used as part of the matrix-vector multiplication
  and the iteration procedure. It is just extra memory.
- `P`: a duck typed matrix to apply the stochastic operator
      the type P must support P*x and P*y    
- `alpha`: the teleportation parameter, alpha must be between 0 and 1.
  This is the probability that the model follows the graph (so the
  problem is hard when alpha is close to 1).
-  `v`: a duck typed vector to apply the personalization
        the type v must support x += v where x is a Vector{T}
        examples of v include a scalar, a sparsevec, or a Vector{T}
- `tol`: the solution tolerance in the error 1-norm.
- `maxiter`: the maximum number of iterations
- `iterfunc`: a function that will be called on each iteration
    
Returns
-------
- `x`: The PageRank vector (just another reference to the input x)

Example
-------
~~~~

~~~~
"""
function pagerank_power!(x::Vector{T}, y::Vector{T},
    P, alpha::T, v, tol::T,
    maxiter::Int, iterfunc::Function) where T
    ialpha = 1.0/(1.0-alpha)
    xinit = x
    _applyv!(x,v,0.,1.0) # iteration number 0
    iterfunc(0,x)
    for iter=1:maxiter
        mul!(y,P,x)
        gamma = 1.0-alpha*KahanSummation.sum_kbn(y)

        delta = 0.
        _applyv!(y,v,alpha,gamma)
        @simd for i=1:length(x)
            @inbounds delta += abs(y[i] - x[i]) # TODO implement Kahan summation
        end
        x,y = y,x # swap
        iterfunc(iter,x)
        if delta*ialpha < tol
            break
        end
    end
    if !(x === xinit) # x is not xinit, so we need to copy it over
        xinit[:] .= x
    end
    xinit # always return xinit
end

"""
`stochastic_mult`
=================

**Expert interface**

This function computes `output` = \$A^T (b .* d)\$, i.e. the matrix \$A\$
tranposed times the vector \$b\$ element-wise multiplied by elements in \$d\$.
The idea is that \$d\$ is the inverse of the row-sums of \$A\$ so that this 
operation computes the (column)-stochastic (or sub-stochastic) product
`output` = \$P b\$ where \$P_{i,j} = A_{j,i} / \$ `degree(j)`.

This function is used by the StochasticMult operators.  
 
Functions
---------
- `stochastic_mult!!(output::StridedVecOrMat, A::SparseMatrixCSC, 
        b::StridedVecOrMat, d::StridedVecOrMat)` 
- `stochastic_mult!!(output::StridedVecOrMat, A::MatrixNetwork, 
        b::StridedVecOrMat, d::StridedVecOrMat)`       

Inputs
------
- `output`: See above 
- `A`: See above
- `b`: See above
- `d`: the vector that will be multipled against b in-place 

Returns
-------
- `output`: the output input is returned as well (as well as updated in place)        
"""
function stochastic_mult end

function stochastic_mult!(output::StridedVecOrMat, A::SparseMatrixCSC, b::StridedVecOrMat, 
    d::StridedVecOrMat)
    
    size(A,1) == size(b,1) || throw(DimensionMismatch())
    size(A,1) == size(output,1) || throw(DimensionMismatch())
    size(output,2) == size(b,2) || throw(DimenensionMismatch())
    nzv = A.nzval
    rv = A.rowval
    fill!(output, zero(eltype(output)))
    for col=1:A.n
        for k=1:size(output,2)
            tmp = zero(eltype(output))
            @inbounds for j=A.colptr[col]:(A.colptr[col+1]-1)
                tmp += (nzv[j])*b[rv[j],k]*d[rv[j]]
            end
            output[col,k] += tmp
        end
    end  
    output  
end    

function stochastic_mult!(output::StridedVecOrMat, A::MatrixNetwork, b::StridedVecOrMat, 
    d::StridedVecOrMat)
    
    # NOTE, for matrix networks, we have already transposed
    # the matrix, so we can do the matvec the other way.
    # that is, this is actually the tranpose one :)
    size(A,1) == size(b,1) || throw(DimensionMismatch())
    size(A,1) == size(output,1) || throw(DimensionMismatch())
    size(output,2) == size(b,2) || throw(DimenensionMismatch())
    fill!(output, zero(eltype(output)))
    # TODO fix the variable names here, they describe 
    # the transpose matrix
    nzv = A.vals
    rv = A.ci
    for col=1:A.n
        for k=1:size(output,2)
            xj = b[col,k]*d[col]
            @inbounds for j=A.rp[col]:(A.rp[col+1]-1)
                output[rv[j],k] += nzv[j]*xj
            end
        end
    end  
    output  
end

"""
`SparseMatrixStochasticMult`
============================

**Expert interface** 

This type is the result of creating
an implicit stochastic operator for
a sparse matrix.
"""
mutable struct SparseMatrixStochasticMult
    d::Vector{Float64}
    A::SparseMatrixCSC     
end 

eltype(op::SparseMatrixStochasticMult) = Float64
size(op::SparseMatrixStochasticMult) = size(op.A)
ndims(op::SparseMatrixStochasticMult) = 2

size(op::SparseMatrixStochasticMult, dim::Integer) = size(op.A,dim)
length(op::SparseMatrixStochasticMult) = length(op.A)

*(op::SparseMatrixStochasticMult, b::AbstractVector{S}) where {S} = 
    mul!(Array{promote_type(Float64,S)}(undef, size(op.A,2)), op, b)
function mul!(output, op::SparseMatrixStochasticMult, b)
    stochastic_mult!(output, op.A, b, op.d)
end

"""
`MatrixNetworkStochasticMult`
============================

**Expert interface** 

This type is the result of creating
an implicit stochastic operator for
a matrix network type.
"""
mutable struct MatrixNetworkStochasticMult
    d::Vector{Float64}
    A::MatrixNetwork     
end 

eltype(op::MatrixNetworkStochasticMult) = Float64
size(op::MatrixNetworkStochasticMult) = size(op.A)
ndims(op::MatrixNetworkStochasticMult) = 2

size(op::MatrixNetworkStochasticMult, dim::Integer) = size(op.A,dim)
length(op::MatrixNetworkStochasticMult) = prod(size(op.A))

*(op::MatrixNetworkStochasticMult, b::AbstractVector{S}) where {S} = 
    mul!(Array{promote_type(Float64,S)}(undef, size(op.A,2)), op, b)
function mul!(output, op::MatrixNetworkStochasticMult, b)
    stochastic_mult!(output, op.A, b, op.d)
end

function _create_stochastic_mult(M::SparseMatrixCSC)
    n = checksquare(M)
    d = Array{Float64}(undef, n, 1)
    sum!(d,M) # compute out-degrees
    for i=1:length(d)
        if d[i]>0.
            d[i] = 1.0/d[i]
        elseif d[i] < 0.
            throw(DomainError(d[i])) 
        end
    end
    Z = SparseMatrixStochasticMult(vec(d),M)
    return Z
end

function _create_stochastic_mult(M::MatrixNetwork)
    A = sparse_transpose(M) # this involves no work...
    n = checksquare(A) 
    d = Array{Float64}(undef,1, n)
    sum!(d,A) # compute out-degrees
    for i=1:length(d)
        if d[i]>0.
            d[i] = 1.0/d[i]
        elseif d[i] < 0.
            throw(DomainError(d[i])) 
        end
    end
    MatrixNetworkStochasticMult(vec(d),M)
end

"""
`pagerank`
---------

PageRank is the stationary distribution of a Markov chain defined
as follows. The behavior of the chain at state i is: 
* with probability \$\alpha\$, randomly transition to an out-neighbor
  of the current node (based on a weighted probabilities if
  the graph has non-negative weights).
* with probability \$1-\alpha\$, jump to a random node chosen anywhere
  in the network. 
* if there are no out-neighbors, then jump anywhere in the network 
  with equal probability. 
    
The solution satisfies a linear system of equations. This function
will solve that linear system to a 1-norm error of `tol` where
`tol` is chosen by default to be the machine precision. 

The solution is always a probability distribution.        

Functions
---------
* `x = pagerank(A::SparseMatrixCSC{V,Int},alpha::Float64)`
* `x = pagerank(A::MatrixNetwork{V},alpha::Float64)`
* `x = pagerank(A,alpha,eps::Float64)` specifies solution tolerance too.        

Inputs
------
- `A`: The sparse matrix or matrix network that you wish to use
  to compute PageRank. In the case of a sparse matrix, it must
  have non-negative values and the values will impact the 
  computation as we will compute a stochastic normalization 
  as part of the algorithm.  
- `alpha`: the teleportation parameter given above.
- `tol`: the tolerance, the default choice is machine precision 
  for floating point, which is more than enough for almost all 
  applications. 

Examples
--------
~~~~
pagerank(A,alpha) 
            # return the standard PageRank vector 
            # with uniform teleportation and alpha=0.85 
            # computed to machine precision            
~~~~
"""
function pagerank end

function pagerank(A,alpha::Float64)
    pagerank(A,alpha,eps(Float64))
end

function pagerank(A,alpha::Float64,tol::Float64)
    _personalized_pagerank_validated(A,alpha,1.0/size(A,1),tol)
end

"""
`seeded_pagerank`
=================

PageRank is the stationary distribution of a Markov chain defined
as follows. The behavior of the chain at state i is: 
* with probability \$\alpha\$, randomly transition to an out-neighbor
  of the current node (based on a weighted probabilities if
  the graph has non-negative weights).
* with probability \$1-\alpha\$, jump according to a distribution
  called the teleportation distribution and given by a vector \$v\$.
  (In the standard case, \$v\$ is the uniform distribution over nodes,
  but this is a parameter).
* if there are no out-neighbors, then transition according to
  the teleportation distribution (this is the strongly-personalized
  problem).
  
When \$v\$ is the uniform vector, then we compute the same thing
as the PageRank call itself.    
    
The solution satisfies a linear system of equations. This function
will solve that linear system to a 1-norm error of `tol` where
`tol` is chosen by default to be the machine precision. 

The solution is always a probability distribution.

Functions
---------
* For any input, A can be either `A::SparseMatrixCSC{V,Int}` for some value
  type V or `A::MatrixNetwork{V}`
* `x = seeded_pagerank(A,alpha::Float64,v::Float64)`
* `x = seeded_pagerank(A,alpha::Float64,v::Int)`
* `x = seeded_pagerank(A,alpha::Float64,v::Set{Int})`
* `x = seeded_pagerank(A,alpha::Float64,v::Dict{Int,Float64})`
* `x = seeded_pagerank(A,alpha::Float64,v::SparseMatrixCSC{Float64})`
* `x = seeded_pagerank(A,alpha::Float64,v::SparseVector{Float64})`
* `x = seeded_pagerank(A,alpha::Float64,v::Vector{Float64})`
* `x = seeded_pagerank(A,alpha::Float64,v::Vector{Float64})`
* `x = seeded_pagerank(A,alpha,v,eps::Float64)` specifies solution tolerance too.        

Inputs
------
- `A`: The sparse matrix or matrix network that you wish to use
  to compute PageRank. In the case of a sparse matrix, it must
  have non-negative values and the values will impact the 
  computation as we will compute a stochastic normalization 
  as part of the algorithm.  
- `alpha`: the teleportation parameter given above.
- `v`: the teleportation distribution vector. This can be
  an integer to teleport to a single node, a set to teleport
  (uniformly) to a set of nodes, a dictionary or sparse vector
  to weight the teleportation. Or a general dense vector. 
- `tol`: the tolerance that we solve the linear system to
  (this is an error guarantee)  

Examples
--------
~~~
seeded_pagerank(A,alpha,5) 
            # return the seeded PageRank vector 
            # with teleportation to node 5 
            # computed to machine precision            
~~~
"""
function seeded_pagerank end

seeded_pagerank(A,alpha,v) = personalized_pagerank(A,alpha,v)
seeded_pagerank(A,alpha,v,tol) = personalized_pagerank(A,alpha,v,tol)

function personalized_pagerank end

function personalized_pagerank(A,alpha::Float64,v)
    personalized_pagerank(A,alpha::Float64,v,eps(Float64))
end

function personalized_pagerank(A,alpha::Float64,v::Float64,tol::Float64)
    if abs(v - 1.0/size(A,1)) >= eps(Float64)
        throw(DomainError(-1))
    end
    _personalized_pagerank_validated(A,alpha,v,tol)
end

function personalized_pagerank(A,alpha::Float64,v::Int,tol::Float64)
    vecv = sparsevec([v], [1.], size(A,1))
    _personalized_pagerank_validated(A,alpha,vecv,tol)
end

function personalized_pagerank(A,alpha::Float64,v::Set{Int},tol::Float64)
    if isempty(v) 
        throw(ArgumentError("the teleportation vector cannot be empty"))
    end
    n = size(A,1)
    if length(v) >= n/3 # TODO validate this size
        # use a dense call
        densev = _densevec(collect(v), ones(Float64, length(v)), size(A,1))
        return personalized_pagerank!(A,alpha,densev,tol) # normalized in the next call 
    else
        sparsev = sparsevec(collect(v), ones(Float64, length(v)), size(A,1))
        return personalized_pagerank!(A,alpha,sparsev,tol) # normalized in the next call
    end 
end

function personalized_pagerank(A,alpha::Float64,v::Dict{Int,Float64},tol::Float64)
    n = size(A,1)
    if isempty(v) 
        throw(ArgumentError("the teleportation vector cannot be empty"))
    end
    if length(v) >= n/3 # TODO validate this size
        # use a dense call
        densev = _densevec(v, size(A,1))
        return personalized_pagerank!(A,alpha,densev,tol) # normalized in the next call 
    else
        sparsev = sparsevec(v, size(A,1))
        return personalized_pagerank!(A,alpha,sparsev,tol) # normalized in the next call
    end
end

personalized_pagerank(A,alpha::Float64,v::SparseMatrixCSC{Float64},tol::Float64) = 
    personalized_pagerank!(A,alpha,deepcopy(v),tol)

function personalized_pagerank!(A,alpha::Float64,v::SparseMatrixCSC{Float64},tol::Float64)
    n = size(A,1)
    if isempty(v) 
        throw(ArgumentError("the teleportation vector cannot be empty"))
    end
    if size(v,1) != n
        throw(ArgumentError(@sprintf("as a sparsevector, v must be n-by-1 where n=%i", n)))
    end
    # This function automatically normalizes the values. 
    vals = nonzeros(v)
    valisum = 1.0/KahanSummation.sum_kbn(vals)
    for ind in eachindex(vals)
        if vals[ind] < 0. 
            throw(DomainError(vals[ind]))
        end
        vals[ind] *= valisum # normalize to sum to 1 
    end
    _personalized_pagerank_validated(A,alpha,v,tol)
end

personalized_pagerank(A,alpha::Float64,v::SparseVector{Float64},tol::Float64) = 
    personalized_pagerank!(A,alpha,deepcopy(v),tol)

function personalized_pagerank!(A,alpha::Float64,v::SparseVector{Float64}, tol::Float64)
    n = size(A,1)
    if isempty(v) 
        throw(ArgumentError("the teleportation vector cannot be empty"))
    end
    if size(v,1) != n
        throw(ArgumentError(@sprintf("as a sparsevector, v must be n-by-1 where n=%i", n)))
    end
    # This function automatically normalizes the values. 
    vals = nonzeros(v)
    valisum = 1.0/KahanSummation.sum_kbn(vals)
    for ind in eachindex(vals)
        if vals[ind] < 0. 
            throw(DomainError(vals[ind]))
        end
        vals[ind] *= valisum # normalize to sum to 1 
    end
    _personalized_pagerank_validated(A,alpha,v,tol)
end


personalized_pagerank(A,alpha::Float64,v::Vector{Float64},tol::Float64) = 
    personalized_pagerank!(A,alpha,deepcopy(v),tol)

function personalized_pagerank!(A,alpha::Float64,v::Vector{Float64},tol::Float64)
    n = size(A,1)
    if isempty(v) 
        throw(ArgumentError("the teleportation vector cannot be empty"))
    end
    if size(v,1) != n
        throw(ArgumentError(@sprintf("as a sparsevector, v must be n-by-1 where n=%i", n)))
    end
    # This function automatically normalizes the values.
    valisum = 1.0/KahanSummation.sum_kbn(v) 
    @inbounds for ind in eachindex(v)
        if v[ind] < 0. 
            throw(DomainError(v[ind]))
        end
        v[ind] *= valisum # normalize to sum to 1 
    end
    _personalized_pagerank_validated(A,alpha,v,tol)
end

function _personalized_pagerank_validated(A,alpha::Float64,v,tol::Float64)
    #x = Vector{Float64}(size(A,1)) # these should suffice, do I still get bugs?
    #y = Vector{Float64}(size(A,1))
    x = zeros(size(A,1))
    y = zeros(size(A,1))
    if alpha < 0. || alpha >= 1.
         throw(ArgumentError(@sprintf("alpha must be in [0,1), but alpha=%f",alpha)))
    end
    maxiter = 2*ceil(Int,log(tol)/log(alpha))
    return pagerank_power!(x, y, _create_stochastic_mult(A), 
                alpha, v, tol, maxiter, _noiterfunc)
end

#------------------------------------------------------------------

"""
`stochastic_heat_kernel_series!`
================================
Compute the stochastic heat kernel, the result
of \$x = \\exp(-t(I - P)) s \$ where \$s\$ is a stochastic
seed vector and \$t\$ is a time parameter and
\$P\$ is a column-stochastic or sub-stochastic
matrix. 

This function will return a vector \$x\$ such that
\$ \\| x - x_{\\mbox{exact}} \\|_1 \\le \$ `eps`, and
this function further guarantees 
\$x_{\\mbox{exact}} - x \\ge 0\$.

The algorithm is just a direct evaluation of the power
series represented by the heat kernel operation.  

** Expert interface, see `seeded_stochastic_heat_kernel` 
for more user-friendly interfaces.**

Inputs
------
-`x`: The solution vector 
-`y`: An intermediate vector 
-`z`: Another intermediate vector
-`P`: a variable that can be used with A_mul_B! 
      in order to apply the stochastic or sub-stochastic matrix
      (column)
-`t`: the value of time in the heat kernel \$0 <= t\$
-`eps`: the solution tolerance \$ 0 < eps < 1 \$
-`maxiter`: the maximum number of series terms to use

Returns
-------
-`x`: the result of the stochastic heat kernel evaluation

Example
-------
~~~~
# TODO Show an example of how to use the expert interface
~~~~

"""
function stochastic_heat_kernel_series!(
    x::Vector{T}, y::Vector{T}, z::Vector{T},
    P, t::T, v, eps::T, 
    maxiter::Int) where T
    
    iexpt = exp(-t)
    _applyv!(y,v,0.,1.) # iteration number 0
    # scale by iexpt
    @simd for i=1:length(x)
        @inbounds x[i] = iexpt*y[i]
    end
    
    eps_exp_t = eps*exp(t)
    err = exp(t)-1.
    coeff = 1.
    
    for k=1:maxiter
        mul!(z,P,y)       # compute z = P*y
        coeff = coeff*t/k  
        BLAS.axpy!(iexpt*coeff, z, x)
        y,z = z,y # swap
        
        err = err - coeff
        if err < eps_exp_t
            break
        end
    end
    x
end



"""
`seeded_stochastic_heat_kernel`
===============================

Compute the stochastic heat kernel, the result
of \$x = \\exp(-t(I - P)) s \$ where \$s\$ is a stochastic
seed vector and \$t\$ is a time parameter and
\$P\$ is a column-stochastic or sub-stochastic
matrix. 

This function will return a vector \$x\$ such that
\$ \\| x - x_{\\mbox{exact}} \\|_1 \\le \$ `eps`, and
this function further guarantees 
\$x_{\\mbox{exact}} - x \\ge 0\$.

Functions
---------
* The input A can be either a `A::SparseMatrixCSC{V,Int}` or `A::MatrixNetwork{V}`
  for some value type `V` 
* `seeded_stochastic_heat_kernel(A,t::Float64,s)` provides the stochastic seed `s` in any
  of the following forms:
* `seeded_stochastic_heat_kernel(A,t,s::Float64)`: s must be the number 1/size(A,1)
  and this corresponds to using the all-ones vector. 
* `seeded_stochastic_heat_kernel(A,t,s::Int)` use the vector with just a single non-zero
* `seeded_stochastic_heat_kernel(A,t,s::Set{Int})` 
* `seeded_stochastic_heat_kernel(A,t,s::Dict{Int,Float64})`
* `seeded_stochastic_heat_kernel(A,t,s::SparseMatrixCSC{Float64,Int})`
* `seeded_stochastic_heat_kernel(A,t,s::SparseVector{Float64,Int})`
* `seeded_stochastic_heat_kernel(A,t,s::Vector{Float64})`      
* `seeded_stochastic_heat_kernel(A,t::Float64,s,eps::Float64)`: specifies 
  the solution tolerance (default `eps=eps(Float64)`)

Inputs
------
-`A`: A sparse matrix with non-negative entries or a matrix network object 
-`t`: the value of time in the heat kernel \$0 <= t\$
-`eps`: the solution tolerance \$ 0 < eps < 1 \$
-`s`: the stochastic seed vector, this will be normalized to be stochastic

Returns
-------
-`x`: the result of the stochastic heat kernel evaluation

Examples
--------
~~~

~~~
"""
function seeded_stochastic_heat_kernel end

function seeded_stochastic_heat_kernel(A,t::Float64,s)
    seeded_stochastic_heat_kernel(A,t::Float64,s,eps(Float64))
end

function seeded_stochastic_heat_kernel(A,t::Float64,s::Float64,tol::Float64)
    if abs(s - 1.0/size(A,1)) >= eps(Float64)
        throw(DomainError(-1))
    end
end

function seeded_stochastic_heat_kernel(A,t::Float64,s::Int,tol::Float64)
    vecv = sparsevec([s], [1.], size(A,1))
    _seeded_heat_kernel_validated(A,t,vecv,tol)
end

function seeded_stochastic_heat_kernel(A,t::Float64,s::Set{Int},tol::Float64)
    if isempty(s) 
        throw(ArgumentError("the seed vector cannot be empty"))
    end
    n = size(A,1)
    if length(s) >= n/3 # TODO validate this size
        # use a dense call
        densev = _densevec(collect(s), ones(Float64, length(s)), size(A,1))
        return seeded_stochastic_heat_kernel!(A,t,densev,tol) # normalized in the next call 
    else
        sparsev = sparsevec(collect(s), ones(Float64, length(s)), size(A,1))
        return seeded_stochastic_heat_kernel!(A,t,sparsev,tol) # normalized in the next call
    end 
end

function seeded_stochastic_heat_kernel(A,t::Float64,s::Dict{Int,Float64},tol::Float64)
    n = size(A,1)
    if isempty(s) 
        throw(ArgumentError("the seed vector cannot be empty"))
    end
    if length(s) >= n/3 # TODO validate this size
        # use a dense call
        densev = _densevec(s, size(A,1))
        return seeded_stochastic_heat_kernel!(A,t,densev,tol) # normalized in the next call 
    else
        sparsev = sparsevec(s, size(A,1))
        return seeded_stochastic_heat_kernel!(A,t,sparsev,tol) # normalized in the next call
    end
end

seeded_stochastic_heat_kernel(A,t::Float64,s::SparseMatrixCSC{Float64},tol::Float64) =
    seeded_stochastic_heat_kernel!(A,t,deepcopy(s),tol)

function seeded_stochastic_heat_kernel!(A,t::Float64,s::SparseMatrixCSC{Float64},tol::Float64)
    n = size(A,1)
    if isempty(s) 
        throw(ArgumentError("the seed vector cannot be empty"))
    end
    if size(s,1) != n
        throw(ArgumentError(@sprintf("as a sparsevector, s must be n-by-1 where n=%i", n)))
    end
    s = copy(s)
    # This function automatically normalizes the values. 
    vals = nonzeros(s)
    valisum = 1.0/KahanSummation.sum_kbn(vals)
    for ind in eachindex(vals)
        if vals[ind] < 0. 
            throw(DomainError(vals[ind]))
        end
        vals[ind] *= valisum # normalize to sum to 1 
    end
    _seeded_heat_kernel_validated(A,t,s,tol)
end

seeded_stochastic_heat_kernel(A,t::Float64,s::SparseVector{Float64},tol::Float64) =
    seeded_stochastic_heat_kernel!(A,t,deepcopy(s),tol)
    
function seeded_stochastic_heat_kernel!(A,t::Float64,s::SparseVector{Float64}, tol::Float64)
    n = size(A,1)
    if isempty(s) 
        throw(ArgumentError("the seed vector cannot be empty"))
    end
    if size(s,1) != n
        throw(ArgumentError(@sprintf("as a sparsevector, s must be n-by-1 where n=%i", n)))
    end
    # This function automatically normalizes the values. 
    vals = nonzeros(s)
    valisum = 1.0/KahanSummation.sum_kbn(vals)
    for ind in eachindex(vals)
        if vals[ind] < 0. 
            throw(DomainError(vals[ind]))
        end
        vals[ind] *= valisum # normalize to sum to 1 
    end
    _seeded_heat_kernel_validated(A,t,s,tol)
end
 # sparsevector

seeded_stochastic_heat_kernel(A,t::Float64,s::Vector{Float64},tol::Float64) =
    seeded_stochastic_heat_kernel!(A,t,deepcopy(s),tol)
      
function seeded_stochastic_heat_kernel!(A,t::Float64,s::Vector{Float64},tol::Float64)
    n = size(A,1)
    if isempty(s) 
        throw(ArgumentError("the seed vector cannot be empty"))
    end
    if size(s,1) != n
        throw(ArgumentError(@sprintf("as a vector, s must be n-by-1 where n=%i", n)))
    end
    # This function automatically normalizes the values.
    valisum = 1.0/KahanSummation.sum_kbn(s) 
    @inbounds for ind in eachindex(s)
        if s[ind] < 0. 
            throw(DomainError(s[ind]))
        end
        s[ind] *= valisum # normalize to sum to 1 
    end
    _seeded_heat_kernel_validated(A,t,s,tol)
end

"""
Returns the number of terms of a Taylor series approximation
to exp(-t(I-P)) sufficient to guarantee eps accuracy in
the one norm when P is a stochastic matrix. 

Since this is actually exp(-t) exp(tP) and P is a stochastic
matrix, if we want 
\$ \\| e^{-t} \\exp(tP) - e^{-t} TaylorPoly(e^t) \\|_1 \\le \\eps \$
then it is equivalent to get
\$ \\| \\exp(tP) -  TaylorPoly(e^t) \\|_1 \\le \\eps e^t  \$
and that is what this code computes. 

Returns -1 if maxdeg is insufficient
"""
function _hk_taylor_degree(t::Float64, eps::Float64, maxdeg::Int) 
    eps_exp_t = eps*exp(t)
    err = exp(t)-1.
    coeff = 1.
    k::Int = 0
    while err > eps_exp_t && k < maxdeg
        k += 1
        coeff = coeff*t/k
        err = err - coeff
    end
    if k==maxdeg && err > eps_exp_t
        return -1
    end 
    return max(k,1)
end

"""
This internal function actually allocates the memory and calls the 
heat kernel algorithm. This could change in the future based on the 
type of graph provided. 
"""
function _seeded_heat_kernel_validated(A,t::Float64,s,tol::Float64)
    x = zeros(size(A,1))
    y = zeros(size(A,1))
    z = zeros(size(A,1))
    if t < 0. || !isfinite(t)
         throw(ArgumentError(@sprintf("t must be in [0,infty), but t=%f",t)))
    end
    #maxiter = 2*ceil(Int,log(tol)/log(alpha))
    bigmaxiter = 10*typemax(Int32) # 21 billion... 
    maxiter = _hk_taylor_degree(t, tol, bigmaxiter)
    if maxiter == -1 
        throw(ArgumentError("the value of t=$t is too large and requires over $bigmaxiter iterations"))
    end  
    return stochastic_heat_kernel_series!(x, y, z,
                 _create_stochastic_mult(A), 
                t, s, tol, maxiter)
end

#=

Just old code

import Base.eltype

type StochasticMult{T}
    n::Int
    mul::Function
end

eltype{T}(op::StochasticMult{T}) = T

ndims(op::StochasticMult) = 2

size(op::StochasticMult, dim::Integer) = (dim == 1) ? op.n :
                                            (dim == 2) ? op.n : 1
size(op::StochasticMult) = (op.m, op.n)

length(op::StochasticMult) = op.m*op.n

*(op::StochasticMult, b) = A_mul_B(op, b)

A_mul_B{R,S}(op::StochasticMult{R}, b::AbstractVector{S}) = A_mul_B!(Array{promote_type(R,S)}(undef, op.m), op, b)

A_mul_B!(output, op::StochasticMult, b) = op.mul(output, b)

function _create_stochastic_mult{V,F}(A::MatrixNetwork{V})
    M = sparse_transpose(A)
    d = Array{F,2}(1,size(M,1))
    sum!(d,M) # compute out-degrees
    for i=1:length(d)
        if d[i]>0.
            d[i] = 1./d[i]
        end
    end
    mul = (x) -> M*(d.*x)
    P = StochasticMult{F}(size(M,1),mul)
end

_create_stochastic_mult{V,F}(A::SparseMatrixCSC{V,Int}) =
    _create_stochastic_mult{V,F}(MatrixNetwork(A))
    
function pagerank{V,FA<:AbstractFloat}(A::MatrixNetwork{V}, alpha::FA)
    F = promote_type(V,FA)
    return pagerank_power(_create_stochastic_mult{V,F}(A),
                F(alpha), F(1./size(A,1)), eps(one(F)))
end

type StochasticMult{T}
    n::Int
    mul::Function
end

eltype{T}(op::StochasticMult{T}) = T
ndims(op::StochasticMult) = 2
size(op::StochasticMult, dim::Integer) = (dim == 1) ? op.n :
                                            (dim == 2) ? op.n : 1
size(op::StochasticMult) = (op.m, op.n)
length(op::StochasticMult) = op.m*op.n
*(op::StochasticMult, b) = A_mul_B(op, b)

A_mul_B{R,S}(op::StochasticMult{R}, b::AbstractVector{S}) = 
    A_mul_B!(Array(promote_type(R,S), op.m), op, b)
A_mul_B!(output, op::StochasticMult, b) = op.mul(output, b)

function _create_stochastic_mult(A::MatrixNetwork{V})
    M = sparse_transpose(A)
    d = Array(Float64, 1, size(M,1))
    sum!(d,M) # compute out-degrees
    for i=1:length(d)
        if d[i]>0.
            d[i] = 1.0/d[i]
        end
    end
    mul = (x) -> M*(d.*x)
    P = StochasticMult{F}(size(M,1),mul)
end

function _handle_pagerank_network(A::SparseMatrixCSC)
    n = checksquare(A)
    d = Array(Float64, n, 1)
    sum!(d,A) # compute out-degrees
    for i=1:length(d)
        if d[i]>0.
            d[i] = 1.0/d[i]
        end
    end
    M = A'
    mul = (x) -> M*(d.*x)
    return StochasticMult{F}(size(M,1),mul)
end

function _create_stochastic_mult(A::SparseMatrixCSC)
    n = checksquare(A)
    vals = nonzeros(A)
    for i=1:length(vals)
        
    end
    d = Array(Float64, n, 1) # needs to be n,1 to get the sum right
    
    sum!(d,A) # compute out-degrees
    for i=1:length(d)
        if d[i]>0.
            d[i] = 1.0/d[i]
        end
    end  
end

"""
`pagerank_power`
----------------

Example
"""
function pagerank_power{F}(P, alpha::F, v, tol::F, _iterfunc=_noiterfunc)
    x = Vector{F}(A.n)
    y = Vector{F}(A.n)
    maxiter = 2*ceil(Int,log(tol)/log(alpha))

    return pagerank_power!(x, y, P, v, tol, maxiter, _iterfunc)
end

function pagerank_power{V,FA<:AbstractFloat}(A::MatrixNetwork{V}, alpha::FA)
    F = promote_type(V,FA)
    P = _create_stochastic_mult{V,F}(A)
end

"""
`pagerank`
--------

Compute a PageRank vector of a graph with adjacency matrix A. This method
solves the linear system of equations

    (I - alpha A' D+) x = (1-alpha) v

such that the solution has a prescribed one-norm error. (The default
is machine precision)

Examples
--------
~~~
A = load_matrix_network("tapir")
pagerank(A,0.85)
pagerank(A,0.85;tol=1e-2) # compute a very low-accuracy PageRank vector
v = rand(size(A,1));
pagerank(A,0.85,v/sum(v)) # personalized PageRank
v = sparsevec([5],[1.],size(A,1)) #
pagerank(A,0.85,v) # personalized PageRank with a sparse right hand side
~~~
"""
function pagerank{V,FA<:AbstractFloat}(A::MatrixNetwork{V}, alpha::FA;
            kwargs...)
    pagerank()
end
function pagerank{V,FA<:AbstractFloat}(A::MatrixNetwork{V}, alpha::FA, v;
            kwargs...)
end
function pagerank{V,FA<:AbstractFloat}(A::SparseMatrixCSC{V}, alpha::FA)
    F = promote_type(V,FA)
    return pagerank_power(_create_stochastic_mult{V,F}(A),
                F(alpha), F(1.0/size(A,1)), eps(one(F)))
end

=#

#=

"""
This function computes the strongly personalized PageRank 
vector of a column sub-stochastic matrix P using the
Neumann series expansion. It computes exactly 

min(maxiter, ...) 

pagerank_neumann
    x - the solution vector, initialized to v
    y - an extra vector of memory
    P - a duck typed matrix to apply the stochastic operator
        the type P must support size(P), eltype(P), P*x 
    alpha - the value of alpha in PageRank
    tol - the solution tolerance
    maxiter - the maximum number of iterations
    iterfunc - 
"""
function pagerank_neumann!{T}(x::Vector{T}, y::Vector{T}, 
    P, alpha::T, v, tol::T, 
    maxiter::Int, iterfunc::Function)
    
    iterfunc(0, x)
    
    for iter=1:maxiter
        A_mul_B!(y,P,x)
        delta = 0.
        @simd for i=1:length(y)
            @inbounds y[i] *= alpha
            @inbounds delta += abs(y[i] - x[i])
        end
        x,y = y,x
        iterfunc(iter,x)
    end
    x
end



=#

