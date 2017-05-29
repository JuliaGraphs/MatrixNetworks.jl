@testset "largest_component" begin
    A = load_matrix_network("dfs_example")
    (Acc,p) = largest_component(A)
    @test size(Acc,1) == 5
    
    (Acc,p) = largest_component(A,true)
    @test size(Acc,1) == 6
    
    A = load_matrix_network("cores_example")
    (Acc1,p1) = largest_component(A)
    (Acc2,p2) = largest_component(A,true)
    @test Acc1 == Acc2
end
