using LinearAlgebra

@testset "generators" begin
    @testset "erdos_renyi" begin
        @test_throws DomainError erdos_renyi_undirected(10,11.)
        @test_throws DomainError erdos_renyi_directed(10,11.)

        n = 100
        avgdegs = range(1., stop=2*log(n), length=100)
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

        @test Matrix(sparse_transpose(havel_hakimi_graph(Int[]))) == zeros(Bool,0,0)
        @test Matrix(sparse_transpose(havel_hakimi_graph([0]))) == zeros(Bool,1,1)
        @test Matrix(sparse_transpose(havel_hakimi_graph([1,1]))) == [0 1; 1 0]
        @test Matrix(sparse_transpose(havel_hakimi_graph([2,2,2]))) == [0 1 1; 1 0 1; 1 1 0]
        six_star = [0 1 1 1 1 1; 1 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0]
        @test Matrix(sparse_transpose(havel_hakimi_graph([5,1,1,1,1,1]))) == six_star

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

    @testset "gpa_graph" begin
        @test typeof(gpa_graph(10,0.5,0.3,2)) == MatrixNetwork{Bool}
        @test typeof(gpa_graph(10,0.5,0.3,2,Val{true})) == MatrixNetwork{Bool}
        @test typeof(gpa_edges!(4,.75,.25,[(1,2)],1)) == Vector{Tuple{Int,Int}}
        @test_throws ArgumentError gpa_edges!(4,.75,.25,[(1,1)],1)
        @test typeof(gpa_edges!(4,.75,.25,[(1,1)],1,Val{true})) == Vector{Tuple{Int,Int}}
        @test is_empty(gpa_graph(0,0.0,0.0,0)) == true
        @test all(diag(sparse_transpose(gpa_graph(10, .3,.4, 3))) .== 0.)
        @test is_undirected(gpa_graph(10, .2, .8,4))
        @test_throws ArgumentError gpa_graph(-1,.3,.4,5)
        @test_throws ArgumentError gpa_graph(8,.8,.4,2)
        @test_throws ArgumentError gpa_graph(10,.3,.4,-2)
        @test size(gpa_graph(12,.5,.4,2),1) == 12
        @test size(gpa_graph(21,.1,.1,3),1) == 21
    end

    @testset "roach_graph" begin
        @test_throws ArgumentError roach_graph(-1)
        @test is_empty(roach_graph(0))
        @test is_connected(roach_graph(5))
        @test is_connected(roach_graph(5, Val{true})[1])
        @test maximum(bfs(roach_graph(5, Val{false}),20)[1]) == 11
        @test maximum(bfs(roach_graph(10),40)[1]) == 21
        @test (bfs(roach_graph(10),40)[1])[1] == 20

        G,xy = roach_graph(5, Val{true})
        filt = vec(all(xy .>= 0,dims = 2))
        Asub = sparse_transpose(G)[filt,filt]
        # Asub should be an antennae...
        @test bfs(Asub,1)[1][5] == 4
        filt = vec(all(xy .<= 0,dims = 2))
        Asub = sparse_transpose(G)[filt,filt]
        # Asub should be an line...
        @test bfs(Asub,1)[1][5] == 4
    end

    @testset "lollipop_graph" begin
        @test_throws ArgumentError lollipop_graph(-1)
        @test_throws ArgumentError lollipop_graph(-1,-1)
        @test_throws ArgumentError lollipop_graph(5,-1,Val{true})

        @test is_connected(lollipop_graph(5))
        @test is_connected(lollipop_graph(5, Val{true})[1])
        @test is_connected(lollipop_graph(5, 10))
        @test is_connected(lollipop_graph(10, 5))

        @test maximum(bfs(lollipop_graph(5,10),1)[1]) == 6
        @test maximum(bfs(lollipop_graph(10,5),1)[1]) == 11

        G,xy = lollipop_graph(6, 5, Val{true})
        filt = vec(all(xy .<= 0, dims = 2))
        Asub = sparse_transpose(G)[filt,filt]
        # Asub should be an tail (or a line graph
        @test bfs(Asub,1)[1][6] == 5
        Asub = sparse_transpose(G)[.~filt,.~filt]
        @test nnz(Asub) == 5*4

    end

    @testset "partial_duplication" begin 
        
        A = sprand(100,100,.2)
        A = max.(A,A')
        B = MatrixNetwork(A)
        @test_throws ArgumentError partial_duplication(B,-1,.1)
        @test_throws ArgumentError partial_duplication(B,1,-.1)
        @test_throws ArgumentError partial_duplication(B,1,1.1)
        @test_throws ArgumentError partial_duplication(MatrixNetwork(sprand(20,10,.1)),1,.1)

        #check edge list creation and conversion
        C,dup_vertices = partial_duplication(B,0,1.0)
        @test length(dup_vertices) == 0
        @test is_undirected(C)
        @test sparse(C) == A

        @inferred partial_duplication(B,100,.5)

        #check symmetry is preserved
        C,dup_vertices = partial_duplication(B,100,.5)
        @test is_undirected(C)
        @test length(dup_vertices) == 100
        
        C,_ = partial_duplication(B,100,.1)
        @test is_undirected(C)

        #check extremal proability case 
        C,_ = partial_duplication(B,100,0.0)
        @test norm(sparse(C),2) == norm(A,2) # no new entries added

    end

    @testset "_get_outedges" begin

        A = sprand(100,100,.2)
        A = max.(A,A')
        B = MatrixNetwork(A)

        @test_throws ArgumentError MatrixNetworks._get_outedges(B,-1)
        @test_throws ArgumentError MatrixNetworks._get_outedges(B,101)
        @inferred MatrixNetworks._get_outedges(B,50)

        row_idx = rand(1:100)
        SA_Is, SA_Vs = findnz(A[row_idx,:])
        MN_Is, MN_Vs = MatrixNetworks._get_outedges(B,row_idx)

        @test MN_Is == SA_Is
        @test MN_Vs == SA_Vs

    end
end
