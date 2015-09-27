function dirclustercoeffs_test()
    (A,xy,labels) = load_matrix_network_metadata("celegans")
    (cc, cccyc, ccmid, ccin, ccout, nf) = dirclustercoeffs(A, true, true)
    return true
end