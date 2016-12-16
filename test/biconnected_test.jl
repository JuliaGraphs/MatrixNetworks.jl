using MatrixNetworks
Output = true
A = load_matrix_network("bicc_example_1")
B = MatrixNetwork(A) 
tic();
bcc = biconnected(B)
toc();
bcc_max = rich_output(B,9)
print(bcc_max)
print("\n")
if (Output)
    print(bcc.bcc_edges)
    print("\n")
    print(bcc.Articulation_points)
    print("\n")
end
