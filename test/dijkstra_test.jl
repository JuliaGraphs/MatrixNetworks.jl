@testset "dijkstra" begin
    (A,xy,labels) = load_matrix_network_metadata("airports")
    A = -A; # fix funny encoding of airport data
    lax = 247; rst = 355
    
    (d,pred) = dijkstra(A,lax)
    @test maximum(d) == 540
    @test pred[end] == lax

    (d,pred) = dijkstra(MatrixNetwork(A),lax)
    @test maximum(d) == 540
    @test pred[end] == lax
    
    @test_throws ErrorException dijkstra(-sparse(1.0I,2,2), 1)
end
