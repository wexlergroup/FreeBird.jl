# show methods

"""
    custom_sort(arr::Vector{Int}, period::Int)

Sorts [1,2,3,4,5,6,7,8] into [1,3,5,7,2,4,6,8] for period = 2. Useful for printing triangular lattices.

"""
function custom_sort(arr::Vector{Int}, period::Int)
    sections = length(arr) ÷ period
    sorted_arr = deepcopy(arr)
    for i in 1:sections
        start = (i-1)*period + 1
        stop = i*period
        sub_arr = arr[start:stop]
        odd_nums = filter(x -> x % 2 != 0, sub_arr)
        even_nums = filter(x -> x % 2 == 0, sub_arr)
        sorted_arr[start:stop] = vcat(odd_nums, even_nums)
    end
    return sorted_arr
end

function print_layer(io::IO, lattice::LatticeSystem{G}, boolvec::Vector{Bool}) where G
    supercell_dimensions = lattice.supercell_dimensions
    # set up zigzag indexing for triangular lattice
    index = custom_sort(collect(1:length(boolvec)), lattice.supercell_dimensions[2]*4)
    ind = 1
    for i in 1:supercell_dimensions[1]
        print(io, "      ")
        if G == TriangularLattice
            if isodd(i)
                print(io, " ")
            end
            for j in 1:supercell_dimensions[2]*length(lattice.basis)
                if boolvec[index[ind]]
                    print(io, "● ")
                else
                    print(io, "○ ")
                end
                ind += 1
            end
        elseif G == SquareLattice
            for j in 1:supercell_dimensions[2]*length(lattice.basis)
                if boolvec[(i-1)*supercell_dimensions[2]*length(lattice.basis) + j]
                    print(io, "● ")
                else
                    print(io, "○ ")
                end
            end
        end
        println(io)
    end
end

function print_lattice(io::IO, lattice::LatticeSystem{G}, boolvec::Vector{Bool}) where G
    if G == GenericLattice
        print(io, boolvec)
        return
    end
    supercell_dimensions = lattice.supercell_dimensions
    if supercell_dimensions[3] == 1
        print_layer(io, lattice, boolvec)
    else
        for k in 1:supercell_dimensions[3]
            # slices the boolvec into layers
            l = supercell_dimensions[1]*supercell_dimensions[2]*length(lattice.basis)
            sec = ((k-1)*l+1):(k*l)
            println(io, "      Layer ", k)
            print_layer(io, lattice, boolvec[sec])
        end
    end
end



# function Base.show(io::IO, lattice::Vector{LatticeSystem{G}}) where G
#     println(io, "Vector{LatticeSystem}(", length(lattice), "):")
#     for (ind, l) in enumerate(lattice)
#         println(io, "[", ind, "] LatticeSystem:")
#         println(io, "    lattice_vectors      : ", l.lattice_vectors)
#         println(io, "    positions            : ", l.positions)
#         println(io, "    supercell_dimensions : ", l.supercell_dimensions)
#         if G == SquareLattice
#             println(io,  "    occupations          : "), print_square_lattice(io, l.supercell_dimensions, l.occupations)
#             println(io,  "    adsorptions          : "), print_square_lattice(io, l.supercell_dimensions, l.adsorptions)
#         elseif G == TriangularLattice
#             println(io,  "    occupations          : "), print_triangular_lattice(io, l.supercell_dimensions, l.occupations)
#             println(io,  "    adsorptions          : "), print_triangular_lattice(io, l.supercell_dimensions, l.adsorptions)
#         else
#             println(io, "    none-square none-triangular lattice")
#             println(io, "    occupations          : ", l.occupations)
#             println(io, "    adsorptions          : ", l.adsorptions)
#         end
#     end
# end
    

# function Base.show(io::IO, lattice::LatticeSystem)
#     println(io, "LatticeSystem:")
#     println(io, "    lattice_vectors      : ", lattice.lattice_vectors)
#     println(io, "    positions            : ", lattice.positions)
#     println(io, "    supercell_dimensions : ", lattice.supercell_dimensions)
#     println(io, "    occupations          : ", lattice.occupations)
#     println(io, "    adsorptions          : ", lattice.adsorptions)
#     println(io, "    neighbors            : ")
#     if length(lattice.neighbors) > 10
#         for i in 1:5
#             println(io, "        site ", i, ": ", "nearest = ", lattice.neighbors[i][1], ", next-nearest = ", lattice.neighbors[i][2])
#         end
#         println(io, "        ⋮")
#         for i in length(lattice.neighbors)-4:length(lattice.neighbors)
#             println(io, "        site ", i, ": ", "nearest = ", lattice.neighbors[i][1], ", next-nearest = ", lattice.neighbors[i][2])
#         end
#     else
#         for i in 1:length(lattice.neighbors)
#             println(io, "        site ", i, ": ", "nearest = ", lattice.neighbors[i][1], ", next-nearest = ", lattice.neighbors[i][2])
#         end
#     end
# end

function Base.show(io::IO, lattice::LatticeSystem{G}) where G
    println(io, "LatticeSystem{$G}:")
    println(io, "    lattice_vectors      : ", lattice.lattice_vectors)
    println(io, "    positions            : ", lattice.positions)
    println(io, "    supercell_dimensions : ", lattice.supercell_dimensions)
    println(io, "    occupations          : "), print_lattice(io, lattice, lattice.occupations)
    println(io, "    adsorptions          : "), print_lattice(io, lattice, lattice.adsorptions)
    # println(io, "    neighbors            : ", lattice.neighbors)
end