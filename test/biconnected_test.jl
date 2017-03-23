function biconnected_test()
    n = 10
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = dropzeros!(sparse(A))
    obj = biconnected_components(B)
    i = 0
    for i = 2:n
        if obj.articulation_points[i]==0
            break
        end
    end
    components = biconnected_components(B;art=false, components=false).number
    @test 2*components == length(obj.map) 
    @test i == n
    @test obj.articulation_points[1]==0
    @test obj.articulation_points[n]==0

    A = load_matrix_network("minnesota")
    B = MatrixNetwork(A)
    components = biconnected_components(B;art=false, components=false).number
    @test components == 142

    A = load_matrix_network("biconnected_example")
    B = MatrixNetwork(A) 
    obj = biconnected_components(B)
    components = biconnected_components(B;art=false, components=false).number
    @test components == 5
    @test obj.articulation_points[5] == 1

    A = load_matrix_network("cores_example")
    B = MatrixNetwork(A)
    obj = biconnected_components(B)
    components = biconnected_components(B;art=false, components=false).number
    @test components == 6
    @test obj.articulation_points[20] == 1
    @test obj.articulation_points[2] == 1
    @test obj.articulation_points[11] == 1
    @test obj.articulation_points[12] == 0

    A = empty_graph()
    obj = biconnected_components(A)
    components = biconnected_components(A;art=false, components=false).number
    @test components == 0

    A = load_matrix_network("clique-10")
    B = MatrixNetwork(A)
    obj = biconnected_components(B)
    components = biconnected_components(B;art=false, components=false).number
    @test components == 1
    return true
end
