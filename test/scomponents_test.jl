function scomponents_test()
    A = load_matrix_network("cores_example")
    cc = scomponents(MatrixNetwork(A))
    
    sizes = [15;1;5]
    components = 3
    if cc.sizes != sizes
        error("scomponents failed")
    end
    if cc.map != 21
        error("scomponents failed")
    end
    if cc.number != components
        error("scomponents failed")
    end
    return true
end