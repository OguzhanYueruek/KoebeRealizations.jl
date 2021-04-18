function johnson_incidence(id::Int)
    if id > 92
        return("Out of range")
    else
    sparse_pm_incidence_matrix = polytope.johnson_solid(id).FACETS_THRU_VERTICES
    # First, we turn the sparse incidence matrix into a dense incidence matrix,
    # then we cast its type into a Julia Bit Matrix, which is the type that koebe_realization() eats.
    return iszero.((-1).+(Matrix(common.dense(sparse_pm_incidence_matrix))))
end
end


function test_johnson_toText(dir)
    #Warmup
    koebe_realization(johnson_incidence(1))

    file_path = dir * "/RationalJohnsonSolids.txt"
    io = open(file_path, "w")
    close(io)
    for i in 1:92
        error = false
        incidence_matrix = johnson_incidence(i)
        result = try
            @timed koebe_realization(incidence_matrix)
        catch
            error = true
            io = open(file_path, "a")
            println(io, "Computation failed for J" * string(i) * "\n")
            close(io)
        end
        if !error
            io = open(file_path, "a")
            println(io, "The vertices of the Koebe realization of J$i are:")
            println(io, result.value)
            println(io, "Computed in "* string(result.time) * "seconds. \n")
            close(io)
        end
    end
end

function test_johnson_toConsole()
    #Warmup
    koebe_realization(johnson_incidence(1))

    for i in 1:92
        error = false
        incidence_matrix = johnson_incidence(i)
        result = try
            @timed koebe_realization(incidence_matrix)
        catch
            error = true
            println("Computation failed for J" * string(i) * "\n")
        end
        if !error
            println("The vertices of the Koebe realization of J$i are:")
            println(result.value)
            println("Computed in "* string(result.time) * "seconds. \n")
        end
    end
end
