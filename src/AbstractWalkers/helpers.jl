# show methods
function print_square_lattice(io::IO, supercell_dimensions::Tuple{Int64, Int64, Int64}, boolvec::Vector{Bool})
    if supercell_dimensions[3] == 1
        for i in 1:supercell_dimensions[1]
            if i != 1
                print(io, "                           ")
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
            if k == 1
                println(io, "Layer ", k)
            else
                println(io, "                           Layer ", k)
            end
            for i in 1:supercell_dimensions[1]
                print(io, "                           ")
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


function print_triangular_lattice(io::IO, supercell_dimensions::Tuple{Int64, Int64, Int64}, boolvec::Vector{Bool})
    if supercell_dimensions[3] == 1
        for i in 1:supercell_dimensions[1]
            if i != 1
                print(io, "                           ")
            end
            if isodd(i)
                print(io, " ")
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
            if k == 1
                println(io, "Layer ", k)
            else
                println(io, "                           Layer ", k)
            end
            for i in 1:supercell_dimensions[1]
                print(io, "                           ")
                if isodd(i)
                    print(io, " ")
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



function Base.show(io::IO, lattice::Vector{LatticeSystem})
    println(io, "Vector{LatticeSystem}(", length(lattice), "):")
    for (ind, l) in enumerate(lattice)
        println(io, "[", ind, "] LatticeSystem:")
        println(io, "    lattice_vectors      : ", l.lattice_vectors)
        println(io, "    positions            : ", l.positions)
        println(io, "    supercell_dimensions : ", l.supercell_dimensions)
        # println(io, "    occupations          : ", l.occupations)
        # println(io, "    adsorptions          : ", l.adsorptions)
        # print(io,  "    occupations          : "), print_square_lattice(io, l.supercell_dimensions, l.occupations)
        # print(io,  "    adsorptions          : "), print_square_lattice(io, l.supercell_dimensions, l.adsorptions)
        print(io,  "    occupations          : "), print_triangular_lattice(io, l.supercell_dimensions, l.occupations)
        print(io,  "    adsorptions          : "), print_triangular_lattice(io, l.supercell_dimensions, l.adsorptions)
    end
end
