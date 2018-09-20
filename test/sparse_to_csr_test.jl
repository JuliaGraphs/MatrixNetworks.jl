@testset "sparse_to_csr" begin
    ei = [1;2;3]
    ej = [3;4;4]
    ev = [8;9;10]
    rp,ci,ai,m = sparse_to_csr(ei,ej,ev)
    
    # more tests added here
    for t = 1:100
        A = sprand(100,80,0.01)
        (rp,ci,ai,m)=sparse_to_csr(A)
        ii = zeros(Int64,length(ai))
        j = zeros(Int64,length(ai))
        a = zeros(Float64,length(ai))
        n = length(rp)-1
        nz = 0
        for cr = 1:n
            for ri = rp[cr]:rp[cr+1]-1
                nz=nz+1
                ii[nz]=cr
                j[nz]=ci[ri]
                a[nz]=ai[ri]
            end
        end
        A2 = sparse(ii,j,a,n,m)
        @test A == A2
    end
end
