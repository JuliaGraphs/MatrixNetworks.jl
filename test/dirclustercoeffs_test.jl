@testset "dirclustercoeffs" begin
    (A,xy,labels) = load_matrix_network_metadata("celegans")
    (cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A, true, true)
    @test length(cc) == 202
    (maxval, maxind) = findmax(cc)
    @test maxind == 113
    (cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A)
    @test length(cc) == 202
    (maxval, maxind) = findmax(cc)
    @test maxind == 113
    (cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A, true)
    @test length(cc) == 202
    (maxval, maxind) = findmax(cc)
    @test maxind == 113
    # negative weights
    @test_throws Exception dirclustercoeffs(-speye(5))
end
