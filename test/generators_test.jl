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

function generators_test()

erdos_renyi_test()
chung_lu_test()
return true

end