using MatrixNetworks

function tridiagnol_test(n::Int64)
    O = ones(Int64,n-1)
    Z = zeros(Int64,n)
    A = Tridiagonal(O,Z,O)
    B = sparse(full(A))
    (components,articulation_points,map) = enrich(B)
    for i = 2:n-1
        if articulation_points[i]==0
            break
        end
    end
    if components != length(map) && i!= n-1 && (articulation_points[1]==1 || articulation_points[n]==1)
        error("biconnected components test failed")
    end
end

function sample_graph_test()
    A = load_matrix_network("biconnected_example")
    B = MatrixNetwork(A) 
    (components,articulation_points) = enrich(A)
    if components!=5 && articulation_point[5]!=1
        error("biconnected components test failed")
    end
end
    
function empty_graph_test()    
    A = empty_graph()
    (components,articulation_points,bcc_edges,component_edges) = enrich(A)
    if components!=0
        error("biconnected components test failed")
    end
end

function sample_clique_test()
    A = load_matrix_network("clique-10")
    B = MatrixNetwork(A)
    (components,articulation_points,map) = enrich(A)
    if components!=1
        error("biconnected components test failed")
    end
end

tridiagnol_test(10)
sample_graph_test()
empty_graph_test()
sample_clique_test()
