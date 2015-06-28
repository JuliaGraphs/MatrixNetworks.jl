# searchdir(path,key) = filter(x->contains(x,key), readdir(path))
# searchdir (generic function with 1 method)
# 
# julia> searchdir("/Users/hudanassar/Documents/A_Research/",".jl")
# 
# #returns the datasets available
function matrix_network_datasets()
    datasets_location = joinpath(Pkg.dir("MatrixNetworks"),"data")
    content = readdir(datasets_location)
    smat_files = filter(x->contains(x,".smat"),content)
    return smat_files
end