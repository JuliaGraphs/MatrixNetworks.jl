# - `cc` = scomponents(A)
# - scomponents(A).num
# - scomponents(A).sizes
# - scomponents(A).map
# - strong_components_map(A)     # if you just want the map
# - enrich(scomponents(A)) # produce additional enriched output
# 
# using MAT
# function scomponents_test()
#     file_path = Pkg.dir("MatrixNetworks/data/bfs_example.mat")
#     file = matopen(file_path)
#     A = read(file,"A")
#     close(file)
#     return bfs(MatrixNetwork(A),1)
# end
