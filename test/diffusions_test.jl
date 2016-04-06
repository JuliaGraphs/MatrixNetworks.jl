include("diffusions.jl")

using Base.Test

n = 10
P = speye(n)
x = zeros(n)
y = zeros(n)
v = rand(n)
v = v/sum(v)
tol = 1e-8
maxiter = 1000
iterfunc = _noiterfunc

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

function _normout(P)
    n = size(P,1)
    colsums = sum(P,1)
    (pi,pj,pv) = findnz(P)
    P = sparse(pi,pj,pv./colsums[pj],n,n)
end
@time P= _normout(P)
x = zeros(n)
y = zeros(n)
v = 1./n
tol = 1e-8
@time x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
