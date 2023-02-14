matrix_to_list_of_list(A::SparseMatrixCSC) = matrix_to_list_of_list(eltype(A.rowval),A)

function matrix_to_list_of_list(::Type{S}, A::SparseMatrixCSC) where S 

    neighbors = Vector{Vector{S}}(undef,size(A,1))

    for i = 1:size(A,1)
        neighbors[i] = A.rowval[A.colptr[i]:A.colptr[i+1]-1]
    end 
    return neighbors

end 

matrix_to_list_of_list(A::MatrixNetwork) = matrix_to_list_of_list(eltype(A.ci),A)


function matrix_to_list_of_list(::Type{S}, A::MatrixNetwork) where S 
    
    neighbors = Vector{Vector{S}}(undef,size(A,1))

    for i = 1:size(A,1)
        neighbors[i] = A.ci[A.rp[i]:A.rp[i+1]-1]
    end 
    return neighbors

end 


function list_of_list_to_sparse_matrix(neighbors::Vector{Vector{S}}) where S 

    Is = S[]
    Js = S[] 
    n = length(neighbors)

    for (v_i,v_i_neighbors) in enumerate(neighbors)
        for v_j in v_i_neighbors
            push!(Is,v_i)
            push!(Js,v_j)
        end 
    end 

    return sparse(Is,Js,1,n,n)
end 