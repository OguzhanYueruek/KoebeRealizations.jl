using KoebeRealizations
using Polymake
using Statistics

"""
    johnson_incidence(n::Int)
    Returns the incidence matrix of the n-th Johnson solid Jn.
"""
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

"""
    are_equal_matrices(v::Vector{Matrix{Float64}})
Checks if all matrices in the vector v are the same, up to a tolerance of 1e-6.
Returns a tuple of boolean and an integer, integer denotes the index of the first different matrix(if v consist of the same matrices, it returns 0)
"""
function are_equal_matrices(v::Vector{Matrix{Float64}})
    a = true
    vertex_mat = v[1]
    first_dif_index = 0
    for k in 2:size(v)[1]
        bitmat = abs.(v[k].-vertex_mat).< 1e-6
        if !prod(prod(bitmat,dims=1))
            a = false
            first_dif_index = k
            break
        end
    end
    return (a,first_dif_index)
end

"""
    test_johnson_toText(dir::String, n::Int, k::Int)

Runs a test of koebe_realization() *n* times for each one of the Johnson solids  from J1 up to J*k*, and prints the result on a text file in the given directory *dir*.
The average computing time, as well as the vertices of the koebe realization is printed. If the vertices of the realization is computed different in any of the *n* trials,
then the set of different vertices is also printed.
"""
function test_johnson_toText(dir,n_trials=3,up_to=92)
    #Warmup
    koebe_realization(johnson_incidence(1))
    file_path = dir * "/RationalJohnsonSolids.txt"
    io = open(file_path, "w")
    close(io)

    for i in 1:up_to
        incidence_matrix = johnson_incidence(i)
        results = []
        for trial in 1:n_trials
            error_flag = false
            result = try
                @timed koebe_realization(incidence_matrix)
            catch
                error_flag = true
                io = open(file_path, "a")
                println(io,"Computation failed for J" * string(i) * " in trial $trial \n")
                close(io)
            end
            if !error_flag
                push!(results, result)
            end
        end
        vertices = [k[1] for k in results]
        timings = [k[2] for k in results]
        time_av = mean(timings)
        io = open(file_path, "a")
        succ_test = length(timings)
        println(io,"Vertices of J" * string(i) * " are computed successfully in $succ_test out of $n_trials trials with an average time of $time_av seconds. \n")
        close(io)
        if !are_equal_matrices(vertices)[1]
            io = open(file_path, "a")
            println(io,"Vertices of J" * string(i) * "in trial " * string(are_equal_matrices(vertices)[2]) * " are different! \n")
            println(io,"Expected \n:" * string(vertices[1]) * "\n Received: \n " * string(vertices[are_equal_matrices(vertices)[2]]) * "\n \n")
            println(io, "<<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>\n")
            close(io)
        else
            io = open(file_path, "a")
            println(io,"Vertices of J" * string(i) * " were same in $succ_test out of $n_trials trials, and are given by \n")
            println(io,vertices[1])
            println(io, "\n <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>> \n")
            close(io)
        end
    end
end



"""
    test_johnson_toText(dir::String, n::Int, k::Int)

Runs a test of koebe_realization() n times for each one of the Johnson solids J1 up to Jk, and prints the result on the console.
The average computing time, as well as the vertices of the koebe realization is printed. If the vertices of the realization is computed different in any of the n trials,
then the set of different vertices is also printed.
"""
function test_johnson_toConsole(n_trials=3,up_to=3)
    #Warmup
    koebe_realization(johnson_incidence(1))

    for i in 1:up_to
        incidence_matrix = johnson_incidence(i)
        results = []
        for trial in 1:n_trials
            error_flag = false
            result = try
                @timed koebe_realization(incidence_matrix)
            catch
                error_flag = true
                println("Computation failed for J" * string(i) * " in trial $trial \n")
            end
            if !error_flag
                push!(results, result)
            end
        end
        vertices = [k[1] for k in results]
        timings = [k[2] for k in results]
        time_av = mean(timings)
        succ_test = length(timings)
        println("Vertices of J" * string(i) * " are computed successfully in $succ_test out of $n_trials trials with an average time of $time_av seconds. \n")
        if !are_equal_matrices(vertices)[1]
            println("Vertices of J" * string(i) * "in trial " * string(are_equal_matrices(vertices)[2]) * " are different! \n")
            println("Expected \n:" * string(vertices[1]) * "\n Received: \n " * string(vertices[are_equal_matrices(vertices)[2]]) * "\n \n")
            println( "<<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>\n")
        else
            println("Vertices of J" * string(i) * " were same in $succ_test out of $n_trials trials, and are given by \n")
            println(vertices[1])
            println("\n <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>\n")
        end
    end
end
