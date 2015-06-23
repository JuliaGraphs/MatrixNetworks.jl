include("../src/readSMAT.jl")
function clustercoeffs_test()
    file_path = Pkg.dir("MatrixNetworks/data/clique-10.smat")
    A = readSMAT(file_path)
    cc = clustercoeffs(MatrixNetwork(A))
    return cc;
end