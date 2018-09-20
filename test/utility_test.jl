using SparseArrays
using LinearAlgebra

@testset "utility" begin
    @testset "random_edge" begin
        @test_throws ArgumentError random_edge(empty_graph(0))
        @test_throws ArgumentError random_edge(empty_graph(5))
        A = sparse([1],[2],1.0,5,5)
        @test random_edge(A)[1:2] == (1,2)
        @test random_edge(A)[1:2] == (1,2)
        @test random_edge(A)[1:2] == (1,2)
        @test random_edge(A)[1:2] == (1,2)
        G = lollipop_graph(5,3)
        A = sparse(G)
        # there are 16 edges, to 100*16*log(16) should randomly generate all of them
        M = zeros(8,8)
        Random.seed!(1) # make it determinstic-ish
        ntrials = 100*16*4
        for i=1:ntrials # log16 = 4 in base 2
            M[random_edge(G)[1:2]...] += 1
        end
        @test sum(M.*A) == ntrials          # this means we got all
        # @test all(findall(!iszero,A.*(M .+ A)) .>= 1) # this means we got all entries
        @test all(x -> iszero(x) || x >= 1, A.*(M .+ A))
        @test std(M[findall(iszero,M)]) <= 100

        n = 16
        A = sparse(1:n-1,2:n,1,n,n)
        # there are 15 edges, to 100*15*log(15) should randomly generate all of them
        M = zeros(n,n)
        Random.seed!(1) # make it determinstic-ish
        ntrials = 100*16*4 # we just use the same one
        for i=1:ntrials # log16 = 4 in base 2
            M[random_edge(A)[1:2]...] += 1
        end
        @test sum(M.*A) == ntrials          # this means we got all
        # @test all(findall(iszero,A.*(M .+ A)) .>= 1) # this means we got all entries
        @test all(x -> iszero(x) || x >= 1, A.*(M .+ A))
        @test std(M[findall(iszero,M)]) <= 100
    end
    @testset "undirected_edges" begin
        @test map(length, undirected_edges(empty_graph(0))) == (0,0)
        @test map(length, undirected_edges(empty_graph(1))) == (0,0)
        @test map(length, undirected_edges(lollipop_graph(5,1))) == (5,5)
        @test undirected_edges(lollipop_graph(5,1)) == ([1,2,3,4,5],[2,3,4,5,6])
        A = sparse([2],[1],1.0,5,5)
        @test map(length, undirected_edges(MatrixNetwork(A))) == (0,0)
    end
    @testset "directed_edges" begin
        @test map(length, directed_edges(empty_graph(0))) == (0,0)
        @test map(length, directed_edges(empty_graph(1))) == (0,0)
        @test map(length, directed_edges(lollipop_graph(5,1))) == (10,10)
        A = sparse([1],[2],1.0,5,5)
        @test map(length, directed_edges(MatrixNetwork(A))) == (1,1)
        @test directed_edges(MatrixNetwork(A)) == ([1],[2])
    end
end