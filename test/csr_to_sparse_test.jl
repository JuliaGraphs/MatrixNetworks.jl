@testset "csr_to_sparse" begin
    ei = [1;2;3]
    ej = [3;4;4]
    ev = [8;9;10]
    (rp,ci,ai,m) = sparse_to_csr(ei,ej,ev)
    (nzi,nzj,nzv) = csr_to_sparse(rp,ci,ai)
    A = sparse(nzi,nzj,nzv,length(rp)-1,maximum(ci))
    
    # more tests added here
    # clique to sparse test
    rp = collect(1:5:26)
    #ci = vec(reshape(repmat(1:5,5,1)',25,1))
    ci = copy(vec(repeat(1:5,5,1)'))
    ai = ones(Int64,25)
    A = csr_to_sparse_matrix(rp,ci,ai,5,5)
    @test Matrix(A) == ones(5,5)
    
    A = csr_to_sparse_matrix(rp,ci,ai)
    @test Matrix(A) == ones(5,5)
    
    # 100 random trials
    for t = 1:100
        A = sprand(100,80,0.01)
        (rp,ci,ai) = sparse_to_csr(A)
        A2 = csr_to_sparse_matrix(rp,ci,ai,100,80)
        @test A == A2
    end
end
