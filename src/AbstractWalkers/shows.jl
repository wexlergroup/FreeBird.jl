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
    elseif sum(lattice.adsorptions) == false
        println(io, "    adsorptions          : no adsorption")
    else
        println(io, "    adsorptions          : "), print_adsorption(io, lattice, lattice.adsorptions)
    end
end

function Base.show(io::IO, lattice::AtomicLattice)
    println(io, typeof(lattice))
    println(io, "    lattice atom         : ", lattice.lattice_atom)
    println(io, "    ASE surface          : ", lattice.surface)
    println(io, "    positions            : ", length(lattice.lattice_positions[:,1]), " grid points")
    println(io, "    supercell_dimensions : ", lattice.supercell_dimensions)
    println(io, "    periodicity          : ", lattice.periodicity)
    println(io, "    adsorbate_atoms      : ", lattice.adsorbate_atoms)
    println(io, "    coverage             : ", coverage(lattice))
    println(io, "    occupied sites       : ", sum(sum, lattice.components), " / ", num_sites(lattice))
    println(io, "    component counts     : ", sum.(lattice.components))
    println(io, "    # nn                 : ", lattice.num_nearest_neighbors)
    println()
end


"""
    view_structure(lattice::AtomicLattice)

Open the lattice in ASE's viewer.

This was `view(lattice::AtomicLattice)`, which defined a function named `view`
inside `AbstractWalkers` and so shadowed `Base.view` for the whole module. It
loaded only because nothing in the module happened to call `Base.view` before
this line — a later `view(A, 1:3)` anywhere in `AbstractWalkers` would have
resolved here and thrown a `MethodError` about `AtomicLattice`.

`view_structure` already exists for `AbstractSystem` and `AtomWalker`
(`helpers.jl`) and is already exported, so this is a method on the right
function rather than a new name. Note it syncs first: `ase_lattice` is a cache
of `occupations` and may be stale.
"""
function view_structure(lattice::AtomicLattice)
    return ase.visualize.view(sync_ase_lattice!(lattice).ase_lattice)
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


"""
    print_lattice_header(io::IO, lattice::AbstractLattice)

Print the geometry summary a live set shows above its walkers.

Dispatched rather than written inline because the two lattice types do not share
a field set: `MLattice` has `lattice_vectors`, `basis` and `cutoff_radii`;
`AtomicLattice` has none of them and carries `lattice_atom`, `type_of_sites` and
a site list instead. The live-set `show` read the `MLattice` names directly,
which was a `FieldError` the moment an `AtomicLattice` could reach it.
"""
function print_lattice_header(io::IO, lattice::MLattice)
    println(io, "    lattice_vectors:      ", lattice.lattice_vectors)
    println(io, "    supercell_dimensions: ", lattice.supercell_dimensions)
    println(io, "    periodicity:          ", lattice.periodicity)
    println(io, "    basis:                ", lattice.basis)
end

function print_lattice_header(io::IO, lattice::AtomicLattice)
    println(io, "    lattice_atom:         ", lattice.lattice_atom)
    println(io, "    adsorbate_atoms:      ", lattice.adsorbate_atoms)
    println(io, "    supercell_dimensions: ", lattice.supercell_dimensions)
    println(io, "    periodicity:          ", lattice.periodicity)
    println(io, "    type_of_sites:        ", lattice.type_of_sites)
    println(io, "    sites:                ", num_sites(lattice))
end

"""
    print_occupation(io::IO, lattice::AtomicLattice)

Print an `AtomicLattice`'s occupancy as a row of 0/1 over `all_sites`.

Needed because `LatticeWalker`'s show path calls `print_occupation` on whatever
configuration it holds, and until `AtomicLattice` could be a walker
configuration at all there was no reason for a method here.
"""
function print_occupation(io::IO, lattice::AtomicLattice)
    merged = zeros(Int, num_sites(lattice))
    for (c, component) in enumerate(lattice.components)
        merged[component] .= c
    end
    print(io, merged)
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
