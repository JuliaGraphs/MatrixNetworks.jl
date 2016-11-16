using MatrixNetworks
Output = false
A = load_matrix_network("bicc_example_1")
B = MatrixNetwork(A) 
bcc = Biconnected(B)
if (Output)
    print(bcc.bcc_edges)
    print("\n")
    print(bcc.Articulation_points)
    print("\n")
end
