function _noiterfunc{T}(iter::Int, x::Vector{T})
end

"""
`_applyv!`
--------

The fully generic function to compute
`x <- alpha*x + gamma*v`

Functions
---------
- `_applyv!(x, v, alpha, gamma)`

"""
function _applyv end

function _applyv!{T}(x::Vector{T}, v, alpha::T, gamma::T)
    x *= alpha
    x += gamma*v
end

function _applyv!{T}(x::Vector{T}, v::T, alpha::T, gamma::T)
    gv = gamma*v
    @simd for i in 1:length(x)
        @inbounds x[i] = alpha*x[i] + gv
    end
end

function _applyv!{T}(x::Vector{T}, v::Vector{T}, alpha::T, gamma::T)
    @simd for i in 1:length(x)
        @inbounds x[i] = alpha*x[i] + gamma*v[i]
    end
end

function _applyv!{T}(x::Vector{T}, v::SparseMatrixCSC{T,Int},
                    alpha::T, gamma::T)
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
This function computes the strongly personalized PageRank
vector of a column sub-stochastic matrix P.

This call allocates no extra memory.

pagerank_power
    x - the solution vector
    y - an extra vector of memory
    P - a duck typed matrix to apply the stochastic operator
        the type P must support P*x
    alpha - the value of alpha in PageRank
    v - a duck typed vector to apply the personalization
        the type v must support x += v where x is a Vector{T}
        examples of v include a scalar, a sparsevec, or a Vector{T}
    tol - the solution tolerance in the error norm.
"""
function pagerank_power!{T}(x::Vector{T}, y::Vector{T},
    P, alpha::T, v, tol::T,
    maxiter::Int, iterfunc::Function)
    ialpha = 1./(1.-alpha)
    xinit = x
    _applyv(x,v,0.,1.) # iteration number 0
    iterfunc(0,x)
    for iter=1:maxiter
        # y = P*x
        A_mul_B!(y,P,x)
        gamma = 1.-alpha*sum_kbn(y)
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
        xinit[:] = x
    end
    xinit # always return xinit
end

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

A_mul_B{R,S}(op::StochasticMult{R}, b::AbstractVector{S}) = A_mul_B!(Array(promote_type(R,S), op.m), op, b)

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

function pagerank_power{F}(P, alpha::F, v, tol::F, _iterfunc=_noiterfunc)
    x = Vector{F}(A.n)
    y = Vector{F}(A.n)
    maxiter = 2*ceil(Int,log(tol)/log(alpha))

    return pagerank_power!(x, y, P, v, tol, maxiter, _iterfunc)
end



function pagerank(A::MatrixNetwork{V}, alpha, v, tol, _iterfunc = _noiterfunc)

end

pagerank_power(P, alpha, v, tol) = pagerank_power(P, alpha, v, tol, _noiterfunc)

function pagerank_power{V,FA<:AbstractFloat}(A::MatrixNetwork{V}, alpha::FA)
    F = promote_type(V,FA)
    P = _create_stochastic_mult{V,F}(A)
end

function pagerank{V,FA<:AbstractFloat}(A::MatrixNetwork{V}, alpha::FA)
    F = promote_type(V,FA)
    return pagerank_power(_create_stochastic_mult{V,F}(A),
                F(alpha), F(1./size(A,1)), eps(one(F)))
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
                F(alpha), F(1./size(A,1)), eps(one(F)))
end


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

