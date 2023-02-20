matrix_to_list_of_list(A::MatrixNetwork) = matrix_to_list_of_list(eltype(A.ci),A)
matrix_to_list_of_list(A::SparseMatrixCSC) = matrix_to_list_of_list(eltype(A.rowval),A)


function matrix_to_list_of_list(::Type{S}, A::MatrixNetwork) where S 
    
    neighbors = Vector{Vector{S}}(undef,A.n)

    for i = 1:A.n
        @inbounds edges,_ = _get_outedges(A,i)
        neighbors[i] = collect(edges)
    end 
    return neighbors

end 

function matrix_to_list_of_list(::Type{S}, A::SparseMatrixCSC) where S 
    
    At = copy(A')

    neighbors = Vector{Vector{S}}(undef,At.n)

    for i = 1:At.n
        @inbounds edges,_ = _get_inedges(At,i)
        neighbors[i] = collect(edges)
    end 
    return neighbors

end 


function list_of_list_to_matrix(::Type{MatrixNetwork{T}}, neighbors::Vector{Vector{S}}) where {S,T}
    return MatrixNetwork(list_of_list_to_matrix(SparseMatrixCSC{T,S},neighbors))
end 

function list_of_list_to_matrix(::Type{SparseMatrixCSC{T1,T2}}, neighbors::Vector{Vector{S}}) where {S,T1,T2} 

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
