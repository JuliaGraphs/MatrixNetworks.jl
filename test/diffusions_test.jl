import KahanSummation
using LinearAlgebra
using SparseArrays

using Test
@testset "diffusions" begin
    function _normout(P::SparseArrays.SparseMatrixCSC{T,Int64}) where T
        n = size(P,1)
        colsums = sum(P,dims=1)
        (pi,pj,pv) = findnz(P)
        P = sparse(pi,pj,pv./colsums[pj],n,n)
    end

    function _normout(P::LinearAlgebra.Adjoint{Float64,SparseArrays.SparseMatrixCSC{T,Int64}}) where T
        n = size(P,1)
        colsums = sum(P,dims=1)
        pj,pi,pv = findnz(P.parent)
        P = sparse(pi,pj,pv./colsums[pj],n,n)
    end

    function prlinsys(A::SparseMatrixCSC{Float64,Int64},alpha::Float64,v::Vector{Float64})
        P = Matrix(_normout(A'))
        n = size(A,1)
        z =  (Matrix(1.0I,n,n) - alpha*P) \ v
        z = z /sum(z)
        return z 
    end

    @testset "pagerank" begin

        A = load_matrix_network("celegans")
        x = pagerank(A, 0.85)
        y = pagerank(A, 0.85, 1e-2)
        
        n = 10
        A = sparse(1.0I,10,10)
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
        
        @test norm(x - sparsevec(collect(1:8),1.0/8,n),1) <= n*eps(Float64)
        

        x = personalized_pagerank(A,0.85,Dict{Int64,Float64}(zip(1:8, ones(8))))
        @test norm(x - sparsevec(collect(1:8),1.0/8,n),1) <= n*eps(Float64)
        
        @test_throws ArgumentError personalized_pagerank(A,0.85,sparsevec(Dict{Int64,Float64}(zip(1:8, ones(8)))))
        
        x = personalized_pagerank(A,0.85,sparsevec(Dict{Int64,Float64}(zip(4:10, ones(7)))))
        @test norm(x - sparsevec(collect(4:10),1.0/7,n),1) <= n*eps(Float64)
        
        x = personalized_pagerank(A, 0.85, sparsevec([5], [2.], 10))
        @test norm(x -xtrue,1) <= n*eps(Float64)
        
        A = sparse(1I,n,n)
        x = personalized_pagerank(A,0.85,5)
        @test norm(x -xtrue,1) <= n*eps(Float64)

        n = 25
        A = LinearAlgebra.fillstored!(copy(sprand(n,n,2/n)), 1)
        x = pagerank(A,0.85)
        xtrue = prlinsys(A,0.85,ones(n))
        @test norm(x -xtrue,1) <= n*eps(Float64)

        x = pagerank(MatrixNetwork(A),0.85)
        @test norm(x -xtrue,1) <= n*eps(Float64)

        n = 10
        A = spdiagm(1=>ones(n-1))
        v = 1.0/n
        tol = 1e-8
        x = pagerank(A,0.85)
        z = zeros(n)
        z[1] = 1.
        z[2] = 1.85
        for ii=3:n
            z[ii] = z[ii-1] + 0.85*(z[ii-1]-z[ii-2])
        end
        z = z/KahanSummation.sum_kbn(z)
        @test norm(x-z,1) <= n*tol
        
    end

    @testset "pagerank_error" begin
        @test_throws DimensionMismatch pagerank(spzeros(6,5),0.85)
        
        @test_throws ArgumentError pagerank(spzeros(5,5),-0.85)
        @test_throws ArgumentError pagerank(spzeros(5,5),1.0)
        @test_throws ArgumentError pagerank(spzeros(5,5),1.01)
        
        @test_throws ArgumentError personalized_pagerank(spzeros(5,5), 0.85, rand(11))
        @test_throws DomainError personalized_pagerank(spzeros(5,5), 0.85, -ones(5))
        @test_throws DomainError personalized_pagerank(spzeros(5,5), 0.85, -sparsevec(ones(5)))    
    end


    @testset "pagerank_alg" begin
        n = 10
        P = sparse(1.0I,n,n)

        x = zeros(n)
        y = zeros(n)
        v = rand(n)
        v = v/sum(v)
        tol = 1e-8
        maxiter = 1000
        iterfunc = MatrixNetworks._noiterfunc

        x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
        @test norm(x-v,1) <= n*tol

        # P = spdiagm(ones(n-1),-1,n,n)
        P = spdiagm(-1=>ones(n-1))
        v .= 0.0
        v[1] = 1.
        x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
        z = zeros(n)
        z[1] = 1.

        for ii=2:n
            z[ii] = 0.85*z[ii-1]
        end
        z = z/KahanSummation.sum_kbn(z)
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

        v = 1.0/n
        tol = 1e-8
        x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,iterfunc)
        z[1] = 1.
        z[2] = 1.85

        for ii=3:n
            z[ii] = z[ii-1] + 0.85*(z[ii-1]-z[ii-2])
        end
        z = z/KahanSummation.sum_kbn(z)

        @test norm(x-z,1) <= n*tol

    end

    @testset "pagerank_perf" begin
        maxiter=500

        n = 1_000
        x = zeros(n)
        y = zeros(n)
        z = zeros(n)
        tol = 1e-8
        v = 1.0/n
        A = LinearAlgebra.fillstored!(copy(sprand(n,n,25/n)),1)

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
        v = 1.0/n
        tol = 1e-8
        A = LinearAlgebra.fillstored!(copy(sprand(n,n,25/n)), 1)
        P = _normout(A')
        
        dt = @elapsed x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
        dt = @elapsed x = pagerank_power!(x,y,P,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
        P2 = MatrixNetworks._create_stochastic_mult(A)
        dt2 = @elapsed y = pagerank_power!(y,z,P2,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
        dt2 = @elapsed y = pagerank_power!(y,z,P2,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)

        #@test dt2 <= 2*dt
        @test norm(x-y) <= n*eps(Float64) 
        
        # now test with matrix networks
        P3 = MatrixNetworks._create_stochastic_mult(MatrixNetwork(A))
        
        y = zeros(n)
        z = zeros(n)
        dt2 = @elapsed y = pagerank_power!(y,z,P3,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
        dt2 = @elapsed y = pagerank_power!(y,z,P3,0.85,v,tol,maxiter,MatrixNetworks._noiterfunc)
        @test norm(x-y) <= n*eps(Float64)
        #@test dt2 <= 2*dt 

        dt = @elapsed x = pagerank_power!(x,y,P,0.85,v,eps(Float64),1000,MatrixNetworks._noiterfunc)
        dt2 = @elapsed y = pagerank(A,0.85)
        @test norm(x-y) <= n*eps(Float64)
        @test dt2 <= 2*dt
    end

    function shkexpm(A,t,v)
        P = Matrix(_normout(A'))
        return exp(-t*(Matrix(1.0I,size(A,1),size(A,1))-P))*v
    end

    @testset "heatkernel" begin

        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,1)
        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,Set([1]))
        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,Dict{Int,Float64}(1 => 1.))
        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,sparsevec([1],[1.],5))
        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,sparse([1],[1],[1.], 5, 1))
        v = zeros(5)
        v[1] = 1.
        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,v)
        seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,1.0/5.0)

        n = 5
        v = rand(5)
        vn = v/sum(v)
        x = seeded_stochastic_heat_kernel(sparse(1.0I,5,5),2.,vn)
        @test norm(x-vn,1) <= n*eps(Float64)
       
        n = 25
        A = LinearAlgebra.fillstored!(copy(sprand(n,n,2/n)), 1)

        v = ones(n)./n
        x = seeded_stochastic_heat_kernel(A,2.,v)
        xtrue = shkexpm(A,2.,v)
        @test norm(x -xtrue,1) <= n*eps(Float64)
        
        x = seeded_stochastic_heat_kernel(MatrixNetwork(A),2.,v)
        @test norm(x -xtrue,1) <= n*eps(Float64)
        
        n = 10
        # A = spdiagm(ones(n-1),-1,n,n)'
        A = copy(spdiagm(-1=>ones(n-1))')
        v = zeros(n)
        v[1] = 1.
        t = 2.
        x = seeded_stochastic_heat_kernel(A,t,v)
        z = zeros(n)
        z[1] = exp(-t)
        for ii=2:n
            z[ii] = z[ii-1]*t./(ii-1.)
        end
        @test norm(x-z,1) <= n*eps(Float64)
        
        @test abs(seeded_stochastic_heat_kernel(spzeros(1,1),2.,1)[1] - exp(-2.)) <= 10*eps(Float64)
        @test abs(seeded_stochastic_heat_kernel(spzeros(1,1),5.,1)[1] - exp(-5.)) <= 10*eps(Float64)
        @test abs(seeded_stochastic_heat_kernel(sparse(1.0I,1,1),2.,1)[1] - 1.) <= 10*eps(Float64) 
    end
    @testset "StochasticMult" begin
        n = 10

        A = MatrixNetwork(LinearAlgebra.fillstored!(copy(sprand(n,n,2/n)), 1))

        MNSM = MatrixNetworks.MatrixNetworkStochasticMult(rand(n), A)
        @test size(MNSM) == (n,n)
        @test size(MNSM, 1) == n
        @test length(MNSM) == n*n
        @test eltype(MNSM) == Float64
        @test ndims(MNSM) == 2

        A = LinearAlgebra.fillstored!(copy(sprand(n,n,2/n)), 1)
        SMSM = MatrixNetworks.SparseMatrixStochasticMult(rand(n), A)
        @test size(SMSM) == (n,n)
        @test size(SMSM, 1) == n
        @test length(SMSM) == n*n
        @test eltype(SMSM) == Float64
        @test ndims(SMSM) == 2
    end
end
