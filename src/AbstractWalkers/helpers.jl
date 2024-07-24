# show methods
function print_lattice(io::IO, lattice::LatticeSystem{G}, boolvec::Vector{Bool}) where G
    if G == GenericLattice
        print(io, boolvec)
        return
    end
    supercell_dimensions = lattice.supercell_dimensions
    if supercell_dimensions[3] == 1
        for i in 1:supercell_dimensions[1]
            print(io, "      ")
            if G == TriangularLattice
                if isodd(i)
                    print(io, " ")
                end
            end
            for j in 1:supercell_dimensions[2]
                if boolvec[(i-1)*supercell_dimensions[2] + j]
                    print(io, "● ")
                else
                    print(io, "○ ")
                end
            end
            println(io)
        end
    else
        for k in 1:supercell_dimensions[3]
            println(io, "      Layer ", k)
            for i in 1:supercell_dimensions[1]
                print(io, "      ")
                if G == TriangularLattice
                    if isodd(i)
                        print(io, " ")
                    end
                end
                for j in 1:supercell_dimensions[2]
                    if boolvec[(i-1)*supercell_dimensions[2] + j + (k-1)*supercell_dimensions[1]*supercell_dimensions[2]]
                        print(io, "● ")
                    else
                        print(io, "○ ")
                    end
                end
                println(io)
            end
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