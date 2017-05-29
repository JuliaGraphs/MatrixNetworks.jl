@testset "generators" begin
    @testset "erdos_renyi" begin
        @test_throws DomainError erdos_renyi_undirected(10,11.)
        @test_throws DomainError erdos_renyi_directed(10,11.)
        
        n = 100
        avgdegs = linspace(1.,2*log(n),100) 
        compsizes = map( (dbar) -> 
                maximum(scomponents(erdos_renyi_undirected(n,dbar)).sizes),
            avgdegs )
            
        @test is_undirected(erdos_renyi_undirected(10,2.))
        
        erdos_renyi_undirected(0, 0.)
        erdos_renyi_undirected(0, 0)
        erdos_renyi_undirected(5, 0)
        erdos_renyi_undirected(5, 0.)
        erdos_renyi_undirected(5, 1.)
        erdos_renyi_undirected(5, 1)
        @test all(diag(sparse_transpose(erdos_renyi_undirected(10, 0.5))) .== 0.)
        erdos_renyi_directed(0, 0)
        erdos_renyi_directed(0, 0)
        erdos_renyi_directed(5, 0.)
        erdos_renyi_directed(5, 1.)
        @test all(diag(sparse_transpose(erdos_renyi_directed(10, 0.5))) .== 0.)
    end

    @testset "chung_lu" begin
        @test_throws ArgumentError chung_lu_undirected([1,2],2)
        @test_throws ArgumentError chung_lu_undirected([1,3])
        @test size(chung_lu_undirected([1,2]),1) == 2
    end

    @testset "havel_hakimi" begin
        @test is_graphical_sequence([1,1]) == true
        @test is_graphical_sequence([0,1,1]) == true
        @test is_graphical_sequence([2,2,2]) == true
        @test is_graphical_sequence([1,1,2]) == true
        @test is_graphical_sequence([1,2,2]) == false
        @test is_graphical_sequence([1,1,3]) == false
        @test is_graphical_sequence([0]) == true
        @test is_graphical_sequence(Int[]) == true
        @test is_graphical_sequence([1,1,1,1]) == true
        
        @test full(sparse_transpose(havel_hakimi_graph(Int[]))) == zeros(Bool,0,0)
        @test full(sparse_transpose(havel_hakimi_graph([0]))) == zeros(Bool,1,1)
        @test full(sparse_transpose(havel_hakimi_graph([1,1]))) == [0 1; 1 0]
        @test full(sparse_transpose(havel_hakimi_graph([2,2,2]))) == [0 1 1; 1 0 1; 1 1 0]
        six_star = [0 1 1 1 1 1; 1 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0] 
        @test full(sparse_transpose(havel_hakimi_graph([5,1,1,1,1,1]))) == six_star
        
        @test_throws ArgumentError is_graphical_sequence([-1])
        @test_throws ArgumentError havel_hakimi_graph([1,2,2])
    end

    @testset "pa_graph" begin
        @test typeof(pa_graph(10,5,5)) == MatrixNetwork{Bool}
        @test typeof(pa_edges!(2,1,[(1,1)])) == Vector{Tuple{Int,Int}}
        @test_throws ArgumentError pa_edges!(5,2,Vector{Tuple{Int,Int}}())
        @test is_empty(pa_graph(0,0,0)) == true
        @test all(diag(sparse_transpose(pa_graph(10, 2, 3))) .== 0.)
        @test is_undirected(pa_graph(10, 12, 3))
        @test_throws ArgumentError pa_graph(-1,10,3)
        @test_throws ArgumentError pa_graph(5,10,-3)
        @test maximum(map(first, pa_edges!(5,1,[(1,2),(2,1)],2))) == 7
        @test size(pa_graph(10,5,5),1) == 10
        @test size(pa_graph(10,-5,5),1) == 10
        @test nnz(sparse_transpose(pa_graph(10,-5,5))) == 20  
    end
end
