 [![Build Result](https://travis-ci.org/nassarhuda/MatrixNetworks.jl.svg?branch=master)](https://travis-ci.org/nassarhuda/MatrixNetworks.jl)
 [![codecov.io](http://codecov.io/github/nassarhuda/MatrixNetworks.jl/coverage.svg?branch=master)](http://codecov.io/github/nassarhuda/MatrixNetworks.jl?branch=master)


# MatrixNetworks
This package consists of a collection of network algorithms.
In short, the major difference between MatrixNetworks.jl and packages like LightGraphs.jl or Graphs.jl is the way graphs are treated.

In [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl), graphs are created through Graph() and DiGraph() which are based on the representation of G as G = (V,E). Similar types exist in [Graphs.jl](https://github.com/JuliaLang/Graphs.jl) (EdgeList, AdjacencyList, IncidenceList, Graph) - this is again based on viewing a graph G as a set of nodes and edges. Our viewpoint is different.

MatrixNetworks is based on the philosophy that there should be no distinction between a matrix and a network - thus the name.

For example, `d,dt,p = bfs(A,1)` computes the bfs distance from the node represented by row 1 to all other nodes of the graph with adjacency matrix A. (A can be of type `SparseMatrixCSC` or `MatrixNetwork`). This representation can be easier to work with and handle.

The package provides documentation with sample runs for all functions - viewable through Juila’s REPL. These sample runs come with sample data, which makes it easier for users to get started on `MatrixNetworks`.


## Package Installation:
##### To install package
```
using Pkg
Pkg.add("MatrixNetworks")
using MatrixNetworks
```

##### Example
```
?bfs
?bipartite_matching
```

##### To run test cases:
```
Pkg.test("MatrixNetworks")
```
## Data available:
##### For a full list of all datasets:
```
matrix_network_datasets()
```
##### Loading data example:
```
load_matrix_network("clique-10")
```

## Some examples:
##### largest_component: Return the largest connected component of a graph
Acc is a sparse matrix containing the largest connected piece of a directed graph A
p is a logical vector indicating which vertices in A were chosen
```
A = load_matrix_network("dfs_example")
Acc,p = largest_component(A)
```

##### clustercoeffs: Compute undirected clustering coefficients for a graph
cc is the clustering coefficients
```
A = load_matrix_network("clique-10")
cc = clustercoeffs(MatrixNetwork(A))
```

##### bfs: Compute breadth first search distances starting from a node in a graph
d is a vector containing the distances of all nodes from node u (1 in the example below)
dt is a vector containing the discover times of all the nodes
pred is a vector containing the predecessors of each of the nodes
```
A = load_matrix_network("bfs_example")
d,dt,pred = bfs(A,1)
```

##### scomponents: Compute the strongly connected components of a graph
```
A = load_matrix_network("cores_example")
sc = scomponents(A)
sc.number #number of connected componenets
sc.sizes #sizes of components
sc.map #the mapping of the graph nodes to their respective connected component
strong_components_map(A) # if you just want the map
sc_enrich = enrich(sc) # produce additional enriched output includes:
sc_enrich.reduction_matrix
sc_enrich.transitive_map
sc_enrich.transitive_order
```
Can work on ei,ej:
```
ei = [1;2;3]
ej = [2;4;1]
scomponents(ei,ej)
```

##### bipartite_matching: Return a maximum weight bipartite matching of a graph
```
ei = [1;2;3]
ej = [3;2;4]
BM = bipartite_matching([10;12;13],ei,ej)
BM.weight
BM.cardinality
BM.match
create_sparse(BM) # get the sparse matrix
edge_list(BM) # get the edgelist
edge_indicator(BM,ei,ej) # get edge indicators
```






