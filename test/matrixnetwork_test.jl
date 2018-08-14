using SparseArrays
using LinearAlgebra

@testset "matrixnetwork" begin
    A = load_matrix_network("dfs_example")
    M = MatrixNetwork(A)
    B = sparse(M)
    @test A == B
    C = sparse_transpose(M)
    @test A' == C
    @test size(M) == size(B)
    @test_throws DomainError size(M, 0)
    @test size(M, 4) == 1
    @test ndims(M) == 2

    load_matrix_network_all("minnesota")
    load_matrix_network_all("U3A")

    #@show issym(sparse([0 1; 0 0]'))
    @test is_undirected(MatrixNetwork(sparse([0 1; 1 0]))) == true
    @test is_undirected(MatrixNetwork(sparse([0 1; 0 1]))) == false
    @test is_undirected(MatrixNetwork(load_matrix_network("dfs_example"))) == false
    @test is_undirected(sparse([0 1; 1 0])) == true

    @test is_connected(havel_hakimi_graph([1,1,1,1])) == false
    @test is_connected(havel_hakimi_graph([2,2,2])) == true
    @test is_connected(sparse([0 1 0; 0 0 1; 1 0 0])) == true
    @test is_connected(spzeros(1,1)) == true
    @test is_connected(empty_graph(0)) == false
    @test is_connected(empty_graph(1)) == true
    @test is_connected(empty_graph(2)) == false

    @test is_empty(empty_graph(0)) == true
    @test is_empty(MatrixNetwork(Int[],Int[],0)) == true
    @test is_empty(MatrixNetwork(Int[1],Int[1])) == false
    @test is_empty(erdos_renyi_undirected(0,0)) == true
    @test is_empty(erdos_renyi_undirected(1,0)) == false

    @testset "matvec" begin
        G = erdos_renyi_directed(5,2)
        y = collect(1:5)
        A = sparse(G)
        @test norm(A*y - G*y,2) <= 10*eps(1.0)
        @test norm(G'*y - A'*y,2) <= 10*eps(1.0)
    end
end

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
