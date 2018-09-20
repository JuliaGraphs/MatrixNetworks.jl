using LinearAlgebra
@testset "biconnected" begin
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = MatrixNetwork(dropzeros!(sparse(A)))
    obj = biconnected_components(B)

    for i = 2:n-1
        @test obj.articulation_points[i]!=0
    end
    components = biconnected_components(B;art=false, components=false).number
    number_of_components = biconnected_components!(B, zeros(Bool,0), zeros(Int64,0))
    @test 2*components == 2*number_of_components == length(obj.map)

    @test obj.articulation_points[1]==0
    @test obj.articulation_points[n]==0
    
    B = dropzeros!(sparse(A))
    obj = biconnected_components(B)

    for i = 2:n-1
        @test obj.articulation_points[i]!=0
    end
    components = biconnected_components(B;art=false, components=false).number

    @test obj.articulation_points[1]==0
    @test obj.articulation_points[n]==0
    
    B = sparse_to_csr(dropzeros!(sparse(A)))
    obj = biconnected_components(B...)
    for i = 2:n-1
        @test obj.articulation_points[i]!=0
    end
    components = biconnected_components(B...;art=false, components=false).number

    @test obj.articulation_points[1]==0
    @test obj.articulation_points[n]==0
    
    B = findnz(dropzeros!(sparse(A)))
    obj = biconnected_components(B[1], B[2])

    for i = 2:n-1
        @test obj.articulation_points[i]!=0
    end
    components = biconnected_components(B[1], B[2];art=false, components=false).number

    @test obj.articulation_points[1]==0
    @test obj.articulation_points[n]==0

    A = load_matrix_network("minnesota")
    B = MatrixNetwork(A)
    components = biconnected_components(B;art=false, components=false).number
    number_of_components = biconnected_components!(B, zeros(Bool,0), zeros(Int64,0))
    @test components == number_of_components == 142

    A = load_matrix_network("biconnected_example")
    B = MatrixNetwork(A) 
    obj = biconnected_components(B)
    components = biconnected_components(B;art=false, components=false).number
    number_of_components = biconnected_components!(B, zeros(Bool,0), zeros(Int64,0))
    @test components == number_of_components == 5
    @test obj.articulation_points[5] == 1

    A = load_matrix_network("cores_example")
    B = MatrixNetwork(A)
    obj = biconnected_components(B)
    components = biconnected_components(B;art=false, components=false).number
    number_of_components = biconnected_components!(B, zeros(Bool,0), zeros(Int64,0))
    @test components == number_of_components == 6
    @test obj.articulation_points[20] == 1
    @test obj.articulation_points[2] == 1
    @test obj.articulation_points[11] == 1
    @test obj.articulation_points[12] == 0

    A = empty_graph()
    obj = biconnected_components(A)
    components = biconnected_components(A;art=false, components=false).number
    number_of_components = biconnected_components!(A, zeros(Bool,0), zeros(Int64,0))
    @test components == number_of_components == 0

    A = load_matrix_network("clique-10")
    B = MatrixNetwork(A)
    obj = biconnected_components(B)
    components = biconnected_components(B;art=false, components=false).number
    number_of_components = biconnected_components!(B, zeros(Bool,0), zeros(Int64,0))
    @test components == number_of_components == 1
end
