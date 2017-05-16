@testset "dirclustercoeffs" begin
    (A,xy,labels) = load_matrix_network_metadata("celegans")
    (cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A, true, true)
    @test length(cc) == 202
    (maxval, maxind) = findmax(cc)
    @test maxind == 113
end
