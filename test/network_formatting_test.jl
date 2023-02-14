using MatrixNetworks: matrix_to_list_of_list, list_of_list_to_sparse_matrix 

@testset "network_formatting" begin

    n = 50 
    
    sp_A = sprand(n,n,.1)
    sp_A = max.(sp_A,sp_A')
    #sp_A.nzval .= 1
    mn_A = MatrixNetwork(sp_A)

    @testset "list of list formatting" begin 


        converted_sp_A = list_of_list_to_sparse_matrix(matrix_to_list_of_list(sp_A))
        @test converted_sp_A.rowval == sp_A.rowval && converted_sp_A.colptr == sp_A.colptr
                #NOTE: using `converted_sp_A == sp_A' fails bc nzval types are different.

        converted_mn_A = MatrixNetwork(list_of_list_to_sparse_matrix(matrix_to_list_of_list(mn_A)))
        @test converted_mn_A.ci == mn_A.ci && converted_mn_A.rp == mn_A.rp 

    end 

end 