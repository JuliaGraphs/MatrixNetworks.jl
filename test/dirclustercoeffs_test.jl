function dirclustercoeffs_test()
    (A,xy,labels) = load_matrix_network_metadata("celegans")
    (cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A, true, true)
    if length(cc) != 202
        error("dirclustercoeffs failed")
    end
    (maxval, maxind) = findmax(cc)
    if maxind != 113
        error("dirclustercoeffs failed")
    end
    
    return true
end