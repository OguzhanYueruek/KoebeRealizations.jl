# This is the part that forms the array of names for Johnson solids
io = open("/home/oguzhan/Documents/Ongoing/KoebeRealizations/KoebeRealizations.jl/src/Solids.txt", "r")
# Read sttring from the text that contains the portion of the polymake docu.
longstr = read(io,String)
split_list = split(longstr,".\n#")
solids = []
for i in split_list
    data = split(i,":")
    if length(data) != 1
        solid_identifier = strip(match(r"J\d+",data[2]).match)
        solid_name = strip(match(r"'\w+'",data[1]).match[2:end-1])
        push!(solids,(solid_identifier,solid_name))
    end
end


name_array = [i[2] for i in solids]


io2 = open("/home/oguzhan/Documents/Ongoing/KoebeRealizations/KoebeRealizations.jl/src/KoebeJohnsonSolids.txt", "w")
close(io2)
for i in 86:length(solids)
    io2 = open("/home/oguzhan/Documents/Ongoing/KoebeRealizations/KoebeRealizations.jl/src/KoebeJohnsonSolids.txt", "a")
    sparse_pm_incidence_matrix = polytope.johnson_solid(solids[i][2]).FACETS_THRU_VERTICES
    write(io2,"The following gives a set of vertices for the Koebe realization of " , solids[i][1], ":\n")

    # First, we turn the sparse incidence matrix into a dense incidence matrix,
    # then we cast its type into a Julia Bit Matrix, which is the type that koebe_realization() eats.
    dense_bit_incidence_matrix = iszero.((-1).+(Matrix(common.dense(sparse_pm_incidence_matrix))))
    print(io2, koebe_realization(transpose(dense_bit_incidence_matrix)), "\n\n")
    close(io2)
end

#open("/home/oguzhan/Documents/Ongoing/KoebeRealizations/KoebeRealizations.jl/src/KoebeJohnsonSolids.txt","w") do f
             # Make sure we write 64bit integer in little-endian byte order
#             write(f,"The following gives a set of vertices for the Koebe realization of")
#end

#19
