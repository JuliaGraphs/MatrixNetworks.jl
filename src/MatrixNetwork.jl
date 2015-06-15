# create type MatrixNetwork
type MatrixNetwork
    n::Int64 # number of columns/rows
    rp::Vector{Int64} # row pointers
    ci::Vector{Int64} # column indices
    vals::Vector{Float64} # corresponding values
end