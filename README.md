# MatrixNetworks

## To install package:
Pkg.clone("https://github.com/nassarhuda/MatrixNetworks.jl.git")

using MatrixNetworks
## To be able to see documentation make sure package Lexicon is installed:
Pkg.add("Lexicon")

using Lexicon

## Example
? bfs

? bipartite_matching

## To run test cases:
Pkg.test("MatrixNetworks")

## Some examples:
### clustercoeffs
```
file_path = Pkg.dir("MatrixNetworks/data/clique-10.smat")
A = readSMAT(file_path)
cc = clustercoeffs(MatrixNetwork(A))
```

### bipartite_matching:
```
ei = [1;2;3]
ej = [3;2;4]
M_out = bipartite_matching([10;12;13],ei,ej)
M_out.weight
M_out.cardinality
M_out.match
MatrixNetworks.create_sparse(bipartite_matching(W)) # get the sparse matrix
MatrixNetworks.edge_list(bipartite_matching(W)) # get the edgelist
MatrixNetworks.edge_indicator(M_out,ei,ej)
```

### scomponents:
```
file_path = Pkg.dir("MatrixNetworks/data/cores_example.smat")
A = readSMAT(file_path)
cc = scomponents(A)
cc.number
cc.sizes
cc.map
strong_components_map(A)     # if you just want the map
enrich(cc) # produce additional enriched output
```

Can work on ei,ej:\n
```
ei = [1;2;3]
ej = [2;4;1]
scomponents(ei,ej)
```

