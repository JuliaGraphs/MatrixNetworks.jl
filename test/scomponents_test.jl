@testset "scomponents" begin
    A = load_matrix_network("cores_example")
    sizes = [15;1;5]
    components = 3
    @testset "MatrixNetwork" begin
        cc = scomponents(MatrixNetwork(A))
        @test cc.sizes == sizes
        @test length(cc.map) == 21
        @test cc.number == components
    end
    @testset "SparseMatrixCSC" begin
        cc = scomponents(A)
        @test cc.sizes == sizes
        @test length(cc.map) == 21
        @test cc.number == components
    end
    @testset "CSR" begin
        cc = scomponents(sparse_to_csr(A)...)
        @test cc.sizes == sizes
        @test length(cc.map) == 21
        @test cc.number == components
    end
    @testset "triplet" begin
        cc = scomponents(findnz(A)[1], findnz(A)[2])
        @test cc.sizes == sizes
        @test length(cc.map) == 21
        @test cc.number == components
    end
    @testset "strong_components_map" begin 
        sci = strong_components_map(MatrixNetwork(A))
        @test strong_components_map(A) == sci
        @test strong_components_map(sparse_to_csr(A)...) == sci
        @test strong_components_map(findnz(A)[1], findnz(A)[2]) == sci
    end
    @testset "empty" begin
        @test scomponents(empty_graph(0)).sizes == []
    end
end
