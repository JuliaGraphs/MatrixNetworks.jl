function erdos_renyi_test()
    @test_throws DomainError erdos_renyi_undirected(10,11.)
    @test_throws DomainError erdos_renyi_directed(10,11.)
    
    n = 100
    avgdegs = linspace(1.,2*log(n),100) 
    compsizes = map( (dbar) -> 
            maximum(scomponents(erdos_renyi_undirected(n,dbar)).sizes),
        avgdegs )
        
    @test is_undirected(erdos_renyi_undirected(10,2.))
    
end

function chung_lu_test() 
    @test_throws ArgumentError chung_lu_undirected([1,2],2)
    @test_throws ArgumentError chung_lu_undirected([1,3])
    @test size(chung_lu_undirected([1,2]),1) == 2
end

function havel_hakimi_test()
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

function generators_test()

erdos_renyi_test()
chung_lu_test()
havel_hakimi_test()
return true

end