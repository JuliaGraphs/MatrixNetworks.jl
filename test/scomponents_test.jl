include("../src/manage_data.jl")
function scomponents_test()
    file_path = Pkg.dir("MatrixNetworks/data/cores_example.smat")
    A = readSMAT(file_path)
    return scomponents(MatrixNetwork(A))
#     scomponents(A).num
#     scomponents(A).sizes
#     scomponents(A).map
#     strong_components_map(A)     # if you just want the map
#     enrich(scomponents(A)) # produce additional enriched output
end