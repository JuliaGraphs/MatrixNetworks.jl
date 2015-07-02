include("../src/manage_data.jl")
function corenums_test()
    file_path = Pkg.dir("MatrixNetworks/data/cores_example.smat")
    A = readSMAT(file_path)
    (d,rt) = corenums(MatrixNetwork(A))
    return (d,rt)
end