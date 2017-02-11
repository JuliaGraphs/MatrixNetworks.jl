using MatrixNetworks
Output = false
n = 10
O = ones(Int64,n-1)
Z = zeros(Int64,n)
A = Tridiagonal(O,Z,O)
B = sparse(full(A))
(cn,articulation_points) = enrich_biconnected(B)
if cn!=n-1 && (articulation_points[1]==1 || articulation_points[n]==1)
    error("biconnected components test failed")
end

A = load_matrix_network("biconnected_example")
B = MatrixNetwork(A) 
(cn,articulation_points) = enrich_biconnected(A)
if cn!=5 && articulation_point[5]!=1
    error("biconnected components test failed")
end

A = empty_graph()
(cn,articulation_points,bcc_edges,component_edges) = enrich_biconnected(A)
if cn!=0
    error("biconnected components test failed")
end
