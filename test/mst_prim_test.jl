include("../src/manage_data.jl")
function mst_prim_test()
    file_path = Pkg.dir("MatrixNetworks/data/airports.smat")
    A = readSMAT(file_path)
    A = -A
    A = max(A,A')
    A = sparse(A)
    T = mst_prim_matrix(A)
    return T
    
    A = load_matrix_network("airports")
A = -A #convert to travel time
A = max(A,A')
A = sparse(A)
(ti,tj,tv) = mst_prim(A)

end