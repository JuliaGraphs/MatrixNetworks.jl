
using MatrixNetworks
#include("../src/diffusions.jl")

using Base.Test

function _normout(P)
    n = size(P,1)
    colsums = sum(P,1)
    (pi,pj,pv) = findnz(P)
    P = sparse(pi,pj,pv./colsums[pj],n,n)
end

function pagerank_test()

    A = load_matrix_network("celegans")
    x = pagerank(A, 0.85)
    #y = pagerank(A, 0.85, 1e.-2)

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