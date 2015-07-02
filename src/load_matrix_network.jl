include("readSMAT.jl")
function load_matrix_network(name::AbstractString)
    basename = joinpath(Pkg.dir("MatrixNetworks"),"data")
    smatfile = joinpath(basename,"$(name).smat")
    if isfile(smatfile)
        return readSMAT(smatfile)
    else
        error(@sprintf "The example datafile '%s' does not seem to exist where it should\n" name)
    end
end

function load_matrix_network_metadata(name::AbstractString)
    basename = joinpath(Pkg.dir("MatrixNetworks"),"data")
    smatfile = joinpath(basename,"$(name).smat")
    meta_xy = joinpath(basename,"$(name).xy")
    meta_labels = joinpath(basename,"$(name).labels")
    if isfile(smatfile)
        if isfile(meta_xy) && isfile(meta_labels)
            xy = readdlm(meta_xy)
            labels = readdlm(meta_labels)
            return (readSMAT(smatfile),xy,labels)
        end
    else
        error(@sprintf "The example datafile '%s' does not seem to exist where it should\n" name)
    end
end