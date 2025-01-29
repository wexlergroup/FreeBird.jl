# show methods

function Base.show(io::IO, walker::AtomWalker{C}) where C
    println(io, "AtomWalker{$C}(")
    println(io, "    configuration      : ", walker.configuration)
    println(io, "    energy             : ", walker.energy)
    println(io, "    iter               : ", walker.iter)
    println(io, "    list_num_par       : ", walker.list_num_par)
    println(io, "    frozen             : ", walker.frozen)
    println(io, "    energy_frozen_part : ", walker.energy_frozen_part,")")
end

function Base.show(io::IO, walker::Vector{AtomWalker{C}}) where C
    println(io, "Vector{AtomWalker{$C}}(", length(walker), "):")
    for (ind, w) in enumerate(walker)
        println(io, "[", ind, "] ", w)
    end
end

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

function Base.show(io::IO, lattice::MLattice{C,G}) where {C,G}
    println(io, typeof(lattice))
    println(io, "    lattice_vectors      : ", lattice.lattice_vectors)
    println(io, "    positions            : ", length(lattice.positions[:,1]), " grid points")
    println(io, "    supercell_dimensions : ", lattice.supercell_dimensions)
    println(io, "    basis                : ", lattice.basis)
    println(io, "    periodicity          : ", lattice.periodicity)
    println(io, "    cutoff radii         : ", length(lattice.cutoff_radii), " nearest neighbors cutoffs ", lattice.cutoff_radii)
    println(io, "    occupations          : ")
    print_occupation(io, lattice)
    if prod(lattice.adsorptions) == true
        println(io, "    adsorptions          : full adsorption")
    else
        println(io, "    adsorptions          : "), print_adsorption(io, lattice, lattice.adsorptions)
    end
end

"""
    merge_components(lattice::MLattice{C}) where C
    
Merges the boolvec of components into a single vector of integers, where each integer represents the component number.
"""
function merge_components(lattice::MLattice{C}) where C
    comp_rep = zeros(Int, prod(lattice.supercell_dimensions)*length(lattice.basis))
    for i in 1:C
        comp_rep += lattice.components[i] * i
    end
    return comp_rep
end

function print_layer(io::IO, lattice::MLattice{C,G}, vec::Union{Vector{Int},Vector{Bool}}) where {C,G}
    if C == 1
        print_layer_single_comp(io, lattice, vec)
    else
        print_layer_multi_comp(io, lattice, vec)
    end
end

function print_layer_multi_comp(io::IO, lattice::MLattice{C,G}, vec::Vector{Int}) where {C,G}
    supercell_dimensions = lattice.supercell_dimensions
    # symlist = [raw"➀", raw"➁", raw"➂", raw"➃", raw"➄", raw"➅", raw"➆", raw"➇", raw"➈", raw"➉"]
    symlist = [raw"➊", raw"➋", raw"➌", raw"➍", raw"➎", raw"➏", raw"➐", raw"➑", raw"➒", raw"➓"]
    # set up zigzag indexing for triangular lattice
    index = custom_sort(collect(1:length(vec)), lattice.supercell_dimensions[2]*4)
    ind = 1
    for i in 1:supercell_dimensions[1]
        print(io, "      ")
        if G == TriangularLattice
            if iseven(i)
                print(io, " ")
            end
            for j in 1:supercell_dimensions[2]*length(lattice.basis)
                if 0 < vec[index[ind]] <= 10
                    print(io, symlist[vec[index[ind]]]*" ")
                elseif vec[index[ind]] > 10
                    print(io, "● ")
                else
                    print(io, "○ ")
                end
                ind += 1
            end
        elseif G == SquareLattice
            for j in 1:supercell_dimensions[2]*length(lattice.basis)
                if 0 < vec[(i-1)*supercell_dimensions[2]*length(lattice.basis) + j] <= 10
                    print(io, symlist[vec[(i-1)*supercell_dimensions[2]*length(lattice.basis) + j]]*" ")
                elseif vec[(i-1)*supercell_dimensions[2]*length(lattice.basis) + j] > 10
                    print(io, "● ")
                else
                    print(io, "○ ")
                end
            end
        end
        println(io)
    end
end

# for single component lattices and adsorptions
function print_layer_single_comp(io::IO, lattice::MLattice{C,G}, boolvec::Vector{Bool}) where {C,G}
    supercell_dimensions = lattice.supercell_dimensions
    # set up zigzag indexing for triangular lattice
    index = custom_sort(collect(1:length(boolvec)), lattice.supercell_dimensions[2]*4)
    ind = 1
    for i in 1:supercell_dimensions[1]
        print(io, "      ")
        if G == TriangularLattice
            if iseven(i)
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


function print_occupation(io::IO, lattice::MLattice{C,G}) where {C,G}
    if G == GenericLattice
        print(io, merge_components(lattice))
        return
    end
    supercell_dimensions = lattice.supercell_dimensions
    if C == 1
        vec = lattice.components[1]
    else
        vec = merge_components(lattice)
    end
    if supercell_dimensions[3] == 1
        print_layer(io, lattice, vec)
    else
        for k in 1:supercell_dimensions[3]
            # slices the boolvec into layers
            l = supercell_dimensions[1]*supercell_dimensions[2]*length(lattice.basis)
            sec = ((k-1)*l+1):(k*l)
            klayer = vec[sec]
            println(io, "     Layer ", k, ":")
            print_layer(io, lattice, klayer)
        end
    end
end


function print_adsorption(io::IO, lattice::MLattice{C,G}, boolvec::Vector{Bool}) where {C,G}
    if G == GenericLattice
        print(io, boolvec)
        return
    end
    supercell_dimensions = lattice.supercell_dimensions
    if supercell_dimensions[3] == 1
        print_layer_single_comp(io, lattice, boolvec)
    else
        for k in 1:supercell_dimensions[3]
            # slices the boolvec into layers
            l = supercell_dimensions[1]*supercell_dimensions[2]*length(lattice.basis)
            sec = ((k-1)*l+1):(k*l)
            println(io, "      Layer ", k)
            print_layer_single_comp(io, lattice, boolvec[sec])
        end
    end
end