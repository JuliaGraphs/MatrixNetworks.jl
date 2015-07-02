function matrix_network_datasets()
    datasets_location = joinpath(Pkg.dir("MatrixNetworks"),"data")
    content = readdir(datasets_location)
    smat_files = filter(x->contains(x,".smat"),content)
    for i = 1:length(smat_files)
        smat_files[i] = smat_files[i][1:end-5]
    end
    return smat_files
end