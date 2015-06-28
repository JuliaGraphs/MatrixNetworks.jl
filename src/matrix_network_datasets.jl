function matrix_network_datasets()
    datasets_location = joinpath(Pkg.dir("MatrixNetworks"),"data")
    content = readdir(datasets_location)
    smat_files = filter(x->contains(x,".smat"),content)
    return smat_files
end