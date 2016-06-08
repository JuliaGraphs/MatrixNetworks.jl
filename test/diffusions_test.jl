
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
    
    x = seeded_pagerank(A,0.85,Set([5]))
    
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
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    x = pagerank(MatrixNetwork(A),0.85)
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    n = 10
    A = spdiagm(ones(n-1),-1,n,n)'
    v = 1./n
    tol = 1e-8
    x = pagerank(A,0.85)
    z = zeros(n)
    z[1] = 1.
    z[2] = 1.85
    for i=3:n
        z[i] = z[i-1] + 0.85*(z[i-1]-z[i-2])
    end
    z = z/sum_kbn(z)
    @test norm(x-z,1) <= n*tol

    
     
    return true 
end

function pagerank_error_test()

    @test_throws DimensionMismatch pagerank(spzeros(6,5),0.85)
    
    @test_throws ArgumentError pagerank(spzeros(5,5),-0.85)
    @test_throws ArgumentError pagerank(spzeros(5,5),1.0)
    @test_throws ArgumentError pagerank(spzeros(5,5),1.01)
    
    @test_throws ArgumentError personalized_pagerank(spzeros(5,5), 0.85, rand(11))
    @test_throws DomainError personalized_pagerank(spzeros(5,5), 0.85, -ones(5))
    @test_throws DomainError personalized_pagerank(spzeros(5,5), 0.85, -sparsevec(ones(5)))    
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

return true
end

function pagerank_perf_test()
    maxiter=500

    n = 1_000
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    tol = 1e-8
    v = 1./n
    A = spones(sprand(n,n,25/n))
    P = _normout(A')
    
    x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
    P2 = MatrixNetworks._create_stochastic_mult(A)
    y = pagerank_power!(y,z,P2,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
    P3 = MatrixNetworks._create_stochastic_mult(MatrixNetwork(A))
    y = pagerank_power!(y,z,P3,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
    
    n = 10_000
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    v = 1./n
    tol = 1e-8
    A = spones(sprand(n,n,25/n))
    P = _normout(A')
    
    dt = @elapsed x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
    P2 = MatrixNetworks._create_stochastic_mult(A)
    dt2 = @elapsed y = pagerank_power!(y,z,P2,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
    
    @test dt2 <= 2*dt
    @test norm(x-y) <= n*eps(Float64) 
    
    # now test with matrix networks
    P3 = MatrixNetworks._create_stochastic_mult(MatrixNetwork(A))
    
    y = zeros(n)
    z = zeros(n)
    dt2 = @elapsed y = pagerank_power!(y,z,P3,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
    
    @test dt2 <= 2*dt
    @test norm(x-y) <= n*eps(Float64)
    
    
    dt = @elapsed x = pagerank_power!(x,y,P,0.85,v,eps(Float64),1000,MatrixNetworks._noiterfunc)
    dt2 = @elapsed y = pagerank(A,0.85)
    @test norm(x-y) <= n*eps(Float64)
    @test dt2 <= 2*dt
     
    
end

function shkexpm(A,t,v)
    P = full(_normout(A'))
    return expm(-t*(eye(size(A,1))-P))*v
end

function heatkernel_test()
    seeded_stochastic_heat_kernel(speye(5),2.,1)
    seeded_stochastic_heat_kernel(speye(5),2.,Set([1]))
    seeded_stochastic_heat_kernel(speye(5),2.,Dict{Int,Float64}(1 => 1.))
    seeded_stochastic_heat_kernel(speye(5),2.,sparsevec([1],[1.],5))
    seeded_stochastic_heat_kernel(speye(5),2.,sparse([1],[1],[1.], 5, 1))
    v = zeros(5)
    v[1] = 1.
    seeded_stochastic_heat_kernel(speye(5),2.,v)
    seeded_stochastic_heat_kernel(speye(5),2.,1./5.)
    
    n = 5
    v = rand(5)
    vn = v/sum(v)
    x = seeded_stochastic_heat_kernel(speye(5),2.,vn)
    @test norm(x-vn,1) <= n*eps(Float64)
   
    n = 25
    A = spones(sprand(n,n,2/n))
    v = ones(n)./n
    x = seeded_stochastic_heat_kernel(A,2.,v)
    xtrue = shkexpm(A,2.,v)
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    x = seeded_stochastic_heat_kernel(MatrixNetwork(A),2.,v)
    @test norm(x -xtrue,1) <= n*eps(Float64)
    
    n = 10
    A = spdiagm(ones(n-1),-1,n,n)'
    v = zeros(n)
    v[1] = 1.
    t = 2.
    x = seeded_stochastic_heat_kernel(A,t,v)
    z = zeros(n)
    z[1] = exp(-t)
    for i=2:n
        z[i] = z[i-1]*t./(i-1.)
    end
    @test norm(x-z,1) <= n*eps(Float64)
    
    @test abs(seeded_stochastic_heat_kernel(spzeros(1,1),2.,1)[1] - exp(-2.)) <= 10*eps(Float64)
    @test abs(seeded_stochastic_heat_kernel(spzeros(1,1),5.,1)[1] - exp(-5.)) <= 10*eps(Float64)
    @test abs(seeded_stochastic_heat_kernel(speye(1,1),2.,1)[1] - 1.) <= 10*eps(Float64) 
end


function diffusions_test()

heatkernel_test()

pagerank_error_test()
pagerank_test()
pagerank_alg_test()
pagerank_perf_test()


return true

end