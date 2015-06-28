include("readSMAT.jl")
function load_matrix_network(name::AbstractString)
    basename = joinpath(Pkg.dir("MatrixNetworks"),"data")
    smatfile = joinpath(basename,"$(name).smat")
    meta_xy = joinpath(basename,"metadata","$(name).xy.smat")
    meta_labels = joinpath(basename,"metadata","$(name).labels.smat")
    if isfile(smatfile)
        if isfile(meta_xy) && isfile(meta_labels)
            xy = readdlm(meta_xy)
            labels = readdlm(meta_labels)
            return (readSMAT(smatfile),xy,labels)
        else
            return readSMAT(smatfile)
        end
    else
        error(@sprintf "The example datafile '%s' does not seem to exist where it should\n" name)
    end
end