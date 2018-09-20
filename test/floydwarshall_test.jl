@testset "floydwarshall" begin
    A = load_matrix_network("all_shortest_paths_example")
    #A.nzval += abs(minimum(A)) + 1 # remove negative edges
    nzvals = nonzeros(A)
    val = abs(minimum(A)) + 1
    for ii=1:length(nzvals)
        nzvals[ii] += val
    end 
    m = size(A,1)
    D2 = zeros(Float64,m,m)
    P2 = zeros(Float64,m,m)
    for ii=1:m
        (d,p) = dijkstra(A,ii)
        D2[ii,:] = d
        P2[ii,:] = p
    end
    
    (D,P) = floydwarshall(A)
    
    @test D == D2
    @test P == P2
    
    (D,P) = floydwarshall(MatrixNetwork(A))
    
    @test D == D2
    @test P == P2
end
