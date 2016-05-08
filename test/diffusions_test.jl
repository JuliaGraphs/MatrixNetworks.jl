
using MatrixNetworks
#include("../src/diffusions.jl")

using Base.Test

function _normout(P)
    n = size(P,1)
    colsums = sum(P,1)
    (pi,pj,pv) = findnz(P)
    P = sparse(pi,pj,pv./colsums[pj],n,n)
end

function prlinsys(A::SparseMatrixCSC{Float64,Int64},alpha::Float64,v::Vector{Float64})
    P = full(_normout(A'))
    n = size(A,1)
    z =  (eye(n) - alpha*P) \ v
    z = z /sum(z)
    return z 
end

function pagerank_test()

    A = load_matrix_network("celegans")
    x = pagerank(A, 0.85)
    y = pagerank(A, 0.85, 1e.-2)
    
    n = 10
    A = speye(10)
    v = rand(n)
    vn = v/sum(v)
    x = personalized_pagerank(A,0.85,v)
    @test norm(x-vn,1) <= n*eps(Float64)
    
    v = 5;
    x = personalized_pagerank(A,0.85,v)
    xtrue = zeros(n)
    xtrue[5] = 1.
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    x = personalized_pagerank(A,0.85,Set([5]))
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    x = personalized_pagerank(A,0.85,Dict{Int64,Float64}(5 => 6.))
    @test norm(x -xtrue,1) <= n*eps(Float64)

    x = personalized_pagerank(A,0.85,Set(collect(1:8)))
    @test norm(x - sparsevec(collect(1:8),1./8,n),1) <= n*eps(Float64)

    x = personalized_pagerank(A,0.85,Dict{Int64,Float64}(zip(1:8, ones(8))))
    @test norm(x - sparsevec(collect(1:8),1./8,n),1) <= n*eps(Float64)
    
    @test_throws ArgumentError personalized_pagerank(A,0.85,sparsevec(Dict{Int64,Float64}(zip(1:8, ones(8)))))
    
    x = personalized_pagerank(A,0.85,sparsevec(Dict{Int64,Float64}(zip(4:10, ones(7)))))
    @test norm(x - sparsevec(collect(4:10),1./7,n),1) <= n*eps(Float64)
    
    x = personalized_pagerank(A, 0.85, sparsevec([5], [2.], 10))
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    A = speye(Int64,n)
    x = personalized_pagerank(A,0.85,5)
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    n = 25
    A = spones(sprand(n,n,2/n))
    x = pagerank(A,0.85)
    xtrue = prlinsys(A,0.85,ones(n))
    @show x
    @show xtrue
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    #A = 
    
     
    return true 
end


function pagerank_alg_test()
n = 10
P = speye(n)
x = zeros(n)
y = zeros(n)
v = rand(n)
v = v/sum(v)
tol = 1e-8
maxiter = 1000
iterfunc = MatrixNetworks._noiterfunc

x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
@test norm(x-v,1) <= n*tol

P = spdiagm(ones(n-1),-1,n,n)
v[:] = 0.
v[1] = 1.
x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
z = zeros(n)
z[1] = 1.
for i=2:n
    z[i] = 0.85*z[i-1]
end
z = z/sum_kbn(z)
@test norm(x-z,1) <= n*tol

tol = 1e-12
x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
@test norm(x-z,1) <= n*tol

tol = 1e-12
x = pagerank_power!(x,y,P,0.85,sparsevec(v),tol,maxiter,iterfunc)
@test norm(x-z,1) <= n*tol


tol = 1e-15
x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
@test norm(x-z,1) <= n*tol

tol = 1e-3
x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,(iter,x) -> @show iter, norm(x,1))
@test norm(x-z,1) <= n*tol

v = 1./n
tol = 1e-8
x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
z[1] = 1.
z[2] = 1.85
for i=3:n
    z[i] = z[i-1] + 0.85*(z[i-1]-z[i-2])
end
z = z/sum_kbn(z)
@test norm(x-z,1) <= n*tol

n = 100_000
P = spones(sprand(n,n,25/n))

x = zeros(n)
y = zeros(n)
v = 1./n
tol = 1e-8
dt = @elapsed x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)

return true
end



function diffusions_test()

pagerank_test()
pagerank_alg_test()


return true

end