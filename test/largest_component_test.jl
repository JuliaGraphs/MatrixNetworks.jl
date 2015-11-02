function largest_component_test()
    A = load_matrix_network("dfs_example")
    (Acc,p) = largest_component(A)
    
    if size(Acc,1) != 5
        error("largest_component failed")
    end
    
    (Acc,p) = largest_component(A,true)
    if size(Acc,1) != 6
        error("largest_component failed")
    end
    
    A = load_matrix_network("cores_example")
    (Acc1,p1) = largest_component(A)
    (Acc2,p2) = largest_component(A,true)
    if !isequal(Acc1,Acc2)
        error("largest_component failed")
    end

    return true
end