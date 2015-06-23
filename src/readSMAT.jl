function readSMAT(filename)
    (rows,header) = readdlm(filename;header=true)
    A = sparse(
               convert(Array{Int64,1},rows[:,1]), 
               convert(Array{Int64,1},rows[:,2]), 
               rows[:,3],
               parseint(Int64,header[1]), 
               parseint(Int64,header[2])
               )
end