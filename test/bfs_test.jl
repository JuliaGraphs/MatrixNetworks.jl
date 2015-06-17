using MAT

file_path = Pkg.dir("MatrixNetworks/data/bfs_example.mat")

file = matopen(file_path)

A = read(file,"A")

close(file)

bfs(MatrixNetwork(A),1)