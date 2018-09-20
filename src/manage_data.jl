using DelimitedFiles

function readSMAT(filename::AbstractString)
    f = open(filename)
    header = readline(f)
    headerparts = split(header)
    nedges = parse(Int,headerparts[3])
    ei = zeros(Int64,nedges)
    ej = zeros(Int64, nedges)
    ev = zeros(Float64, nedges)
    @inbounds for i = 1:nedges
        curline = readline(f)
        parts = split(curline)
        ei[i] = parse(Int, parts[1])+1
        ej[i] = parse(Int, parts[2])+1
        ev[i] = parse(Float64, parts[3])
    end
    close(f)
    A = sparse(ei, ej, ev,
               parse(Int,headerparts[1]), 
               parse(Int,headerparts[2])
               )
    return A
end



mutable struct MatrixNetworkMetadata
    A::SparseMatrixCSC{Int64,Int64}
    labels::Vector{AbstractString}
    xy::Array{Float64,2}
    source::AbstractString
end

function load_matrix_network_all(name::AbstractString)
    A = load_matrix_network(name)
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    
    meta_source = joinpath(pathname,"$(name).source")
    if isfile(meta_source)
        # source = open(readstring, meta_source)
        source = open(s->read(s,String),meta_source)
    else
        source = "(None given)"
    end
    
    meta_xy = joinpath(pathname,"$(name).xy")
    if isfile(meta_xy)
        xy = readdlm(meta_xy)
    else
        xy = zeros(0,2)
    end
    
    meta_labels = joinpath(pathname,"$(name).labels")
    if isfile(meta_labels)
        labels = open(readlines, meta_labels)
    else
        labels = map(x -> @sprintf("%i",x), 1:size(A,1))
    end

    return MatrixNetworkMetadata(A,labels,xy,source)
end


function load_matrix_network(name::AbstractString)
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    smatfile = joinpath(pathname,"$(name).smat")
    if isfile(smatfile)
        return readSMAT(smatfile)
    else
        error(@sprintf "The example datafile '%s' does not seem to exist where it should\n" name)
    end
end

function load_matrix_network_metadata(name::AbstractString)
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    smatfile = joinpath(pathname,"$(name).smat")
    meta_xy = joinpath(pathname,"$(name).xy")
    meta_labels = joinpath(pathname,"$(name).labels")
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

function matrix_network_datasets()
    datasets_location = joinpath(dirname(dirname(@__FILE__)),"data")
    content = readdir(datasets_location)
    cc = map(i->match(r".smat",content[i]),1:length(content))
    smat_files = content[findall(typeof.(cc).!=Nothing)]
    for i = 1:length(smat_files)
        smat_files[i] = smat_files[i][1:end-5]
    end
    return smat_files
end

