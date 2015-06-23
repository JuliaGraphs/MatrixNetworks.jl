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
```
    file_path = Pkg.dir("MatrixNetworks/data/clique-10.smat")
    A = readSMAT(file_path)
    cc = clustercoeffs(MatrixNetwork(A))
```

