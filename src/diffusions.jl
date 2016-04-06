
function _noiterfunc{T}(iter::Int, x::Vector{T})
end

"""
The fully generic function
x <- alpha*x + gamma*v
"""
function _applyv{T}(x::Vector{T}, v, alpha::T, gamma::T)
    x *= alpha
    x += gamma*v
end

function _applyv{T}(x::Vector{T}, v::T, alpha::T, gamma::T)
    gv = gamma*v
    @simd for i in 1:length(x)
        @inbounds x[i] = alpha*x[i] + gv
    end
end

function _applyv{T}(x::Vector{T}, v::Vector{T}, alpha::T, gamma::T)
    @simd for i in 1:length(x)
        @inbounds x[i] = alpha*x[i] + gamma*v[i]
    end
end

function _applyv{T}(x::Vector{T}, v::SparseMatrixCSC{T,Int}, 
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
        _applyv(y,v,alpha,gamma)
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



