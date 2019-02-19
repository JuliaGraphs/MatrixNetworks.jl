@testset "triangles" begin
    A = load_matrix_network("clique-10")
    mytriangles = triangles(A)
    z = collect(mytriangles)
    @test length(z) == 120
    ei,ej,ek = unzip_triangles(z)
    @test length(ei) == 120

    A = erdos_renyi_undirected(100,0.2)
    z = collect(triangles(A;symmetries = true))
    @test sum(Diagonal(sparse(A)^3)) == length(z)
end
