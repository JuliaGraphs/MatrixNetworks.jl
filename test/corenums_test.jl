using MAT
function corenums_test()
    file_path = Pkg.dir("MatrixNetworks/data/cores_example.mat")
    file = matopen(file_path)
    A = read(file,"A")
    close(file)
    (d,rt) = corenums(MatrixNetwork(A))
    return (d,rt)
end