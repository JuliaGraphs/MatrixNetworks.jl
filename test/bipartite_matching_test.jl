using LinearAlgebra

@testset "bipartite_matching_test" begin
    W = sprand(10,8,0.5)
    bipartite_matching(W)
    bipartite_matching([10;12;13],[1;2;3],[3;2;4])
    
    A = [5 1 1 1 1;
     1 5 1 1 1;
     1 1 5 1 1;
     1 1 1 5 1;
     1 1 1 1 5;
     1 1 1 1 1]
    A = sparse(A)
    M1 = bipartite_matching(A)
    (m1,m2) = edge_list(M1)
     
    @test M1.weight == 25/maximum(abs.(A)) && isequal(m1,m2)
    
    A = ones(Int64,6,5) - Matrix(I,6,5)
    A = sparse(A')
    (ei,ej,ev) = csr_to_sparse(A.colptr,A.rowval,A.nzval)
    ai = [1;2;3;4;5;ei]
    aj = [1;2;3;4;5;ej]
    av = [5;5;5;5;5;ev]
    
    M2 = bipartite_matching(av,ai,aj)
    mi = MatrixNetworks.edge_indicator(M2,ai,aj)
    mitrue = zeros(Int64,length(av))
    mitrue[1:5] .= 1
    
    @test M2.weight == 25/maximum(abs.(av)) && isequal(m1,m2) && sum(mi) ==5 && isequal(mi,mitrue)
    
    M2 = bipartite_matching(av,ai,aj,maximum(ai),maximum(aj))
    mi = MatrixNetworks.edge_indicator(M2,ai,aj)
    mitrue = zeros(Int64,length(av))
    mitrue[1:5] .= 1
    
    @test M2.weight == 25/maximum(abs.(av)) && isequal(m1,m2) && sum(mi) ==5 && isequal(mi,mitrue)
end


@testset "bipartite_cardinality_matching_test" begin
    mo_car = bipartite_cardinality_matching([1,2,3],[3,2,4])
    mo_car_sorted = bipartite_cardinality_matching([1,2,3],[3,2,4]; ei_sorted=true)
    @test mo_car.cardinality == 3
    @test mo_car_sorted.cardinality == 3
    @test allunique(mo_car_sorted.match)

    # max cardinality not as big as m
    ei = [1,1,3,3,3,3,1,2,4,4,5,6,6,7,8,8]
    ej = [2,3,2,3,4,5,6,1,4,5,6,5,6,2,7,8]
    w = ones(length(ei))
    mo_car = bipartite_cardinality_matching(ei, ej)
    mo_max = bipartite_matching(w, ei, ej)
    @test mo_car.cardinality == mo_max.cardinality
    @test mo_car.cardinality == 7
    # unmatched vertices have a 0 at the position
    matched = filter(v->v!=0, mo_car.match)
    @test length(matched) == mo_car.cardinality
    @test allunique(matched)

    # max cardinality now same as m
    ei = [1,1,3,3,3,3,1,2,4,4,5,6,6,7,8,8,4]
    ej = [2,3,2,3,4,5,6,1,4,5,6,5,6,2,7,8,8]
    w = ones(length(ei))
    mo_car = bipartite_cardinality_matching(ei, ej)
    mo_max = bipartite_matching(w, ei, ej)
    @test mo_car.cardinality == mo_max.cardinality
    @test mo_car.cardinality == 8
    @test allunique(mo_car.match)

    # m > n
    ei = [1,1,3,3,3,3,1,2,4,4,5,6,6,7,8,8,4,9]
    ej = [2,3,2,3,4,5,6,1,4,5,6,5,6,2,7,8,8,1]
    w = ones(length(ei))
    mo_car = bipartite_cardinality_matching(ei, ej)
    mo_max = bipartite_matching(w, ei, ej)
    @test mo_car.cardinality == mo_max.cardinality
    @test mo_car.cardinality == 8
    # unmatched vertices have a 0 at the position
    matched = filter(v->v!=0, mo_car.match)
    @test length(matched) == mo_car.cardinality
    @test allunique(matched)

    # having a sparse matrix
    ei = [1,1,3,3,3,3,1,2,4,4,5,6,6,7,8,8,4]
    ej = [2,3,2,3,4,5,6,1,4,5,6,5,6,2,7,8,8]
    w = ones(length(ei))
    A = sparse(ei, ej, w)
    mo_car = bipartite_cardinality_matching(A)
    @test mo_car.cardinality == 8
    @test allunique(mo_car.match)
end