# A = load_matrix_network_graph("clique-10")
# (A,meta) = load_matrix_network_data("clique-10") # load additional metadata
# matrix_network_datasets() # Provide a list of datasets
# The implementation of load_matrix_network_graph would just be really simple:
# 
# function load_matrix_network_graph(name::AbstractString)
# 
#    basename = joinpath(dirname(@__FILE__), "..", "data", name)
# 
#     smatname = joinpath(basename, string(dataset_name, ".smat"))
#     if isfile(rdaname)
#         return readSMAT(read_rda(rdaname)[dataset_name])
#     else
#         error(@sprintf "The example datafile '%s' does not seem to exist where it should\n" name)
#     end
# 
# end

# using MatrixNetworks
# using Base.Test
# using Compat
# 
# all_tests = ["dfs_test", "bfs_test", "bipartite_matching_test", "clustercoeffs_test", 
#              "scomponents_test"]
# 
# for t in all_tests
#     test_path = joinpath(Pkg.dir("MatrixNetworks"), "test", "$(t).jl")
#     println("running $(test_path) ...")
#     include(test_path)
# end
# 
# test_path = joinpath(Pkg.dir("MatrixNetworks"), "test", "$(t).jl")
#     println("running $(test_path) ...")
#     include(test_path)
# 
# function readSMAT(filename)
#     (rows,header) = readdlm(filename;header=true)
#     A = sparse(
#                convert(Array{Int64,1},rows[:,1]), 
#                convert(Array{Int64,1},rows[:,2]), 
#                rows[:,3],
#                parse(Int,header[1]), 
#                parse(Int,header[2])
#                )
# end


function load_matrix_network(name::ASCIIString)
    basename = joinpath(Pkg.dir("MatrixNetworks"),"data")
    smatfile = joinpath(basename,"$(name).smat")
    meta_xy = joinpath(basename,"metadata","$(name).xy.smat")
    meta_labels = joinpath(basename,"metadata","$(name).labels.smat")
    if isfile(smatfile)
        if isfile(meta_xy) && isfile(meta_labels)
            xy = readdlm(meta_xy)
            labels = readdlm(meta_lables)
            return (readSMAT(basename),xy,labels)
        else
            return readSMAT(basename)
        end
    else
        error(@sprintf "The example datafile '%s' does not seem to exist where it should\n" name)
    end
end