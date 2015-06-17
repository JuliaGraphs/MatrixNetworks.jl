using MAT

# file_path = Pkg.dir("MatrixNetworks/data/bfs_example.mat")

file = matopen("../data/bfs_example.mat")

A = read(file,"A")

close(file)

bfs(MatrixNetwork(A),1)