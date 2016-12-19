using MatrixNetworks
Output = false
A = load_matrix_network("biconnected_example")
#A = empty_graph(5)
#A = empty_graph()
B = MatrixNetwork(A) 
tic();
bcc = biconnected(B)
toc();
(bcc_edges,component_edges) = enrich_biconnected(B,0)
print(bcc_edges)
print("\n")
print(component_edges)
if (Output)
    print(bcc.biconnected_components_label)
    print("\n")
    print(bcc.articulation_points)
    print("\n")
end
