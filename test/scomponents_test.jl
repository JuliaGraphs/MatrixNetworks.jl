function scomponents_test()
    A = load_matrix_network("cores_example")
    cc = scomponents(MatrixNetwork(A))
    
    sizes = [15;1;5]
    components = 3
    if cc.sizes != sizes
        error("scomponents failed")
    end
    if length(cc.map) != 21
        error("scomponents failed")
    end
    if cc.number != components
        error("scomponents failed")
    end
    
    # test that SparseMatrixCSC works
    cc = scomponents(A)
    
    sizes = [15;1;5]
    components = 3
    if cc.sizes != sizes
        error("scomponents failed")
    end
    if length(cc.map) != 21
        error("scomponents failed")
    end
    if cc.number != components
        error("scomponents failed")
    end
    
    # test that a csr matrix works
    cc = scomponents(sparse_to_csr(A)...)
    
    sizes = [15;1;5]
    components = 3
    if cc.sizes != sizes
        error("scomponents failed")
    end
    if length(cc.map) != 21
        error("scomponents failed")
    end
    if cc.number != components
        error("scomponents failed")
    end
    return true
    
    # test that I,J tuple works
    cc = scomponents(findnz(A)[1], findnz(A)[2])
    
    sizes = [15;1;5]
    components = 3
    if cc.sizes != sizes
        error("scomponents failed")
    end
    if length(cc.map) != 21
        error("scomponents failed")
    end
    if cc.number != components
        error("scomponents failed")
    end

    sci = strong_components_map(MatrixNetwork(A))
    @test strong_components_map(A) == sci
    @test strong_components_map(sparse_to_csr(A)...) == sci
    @test strong_components_map(findnz(A)[1], findnz(A)[2]) == sci
end
