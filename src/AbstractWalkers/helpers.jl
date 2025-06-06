# convenience functions for handling atomistic walkers

"""
    split_components(at::AbstractSystem, list_num_par::Vector{Int})

Split the system into components based on the number of particles in each component.

# Arguments
- `at::AbstractSystem`: The system to split.
- `list_num_par::Vector{Int}`: The number of particles in each component.

# Returns
- `components`: An array of `FastSystem` objects representing the components of the system.

"""
function split_components(at::AbstractSystem, list_num_par::Vector{Int})
    components = Array{FastSystem}(undef, length(list_num_par))
    comp_cut = vcat([0],cumsum(list_num_par))
    comp_split = [comp_cut[i]+1:comp_cut[i+1] for i in 1:length(list_num_par)]
    for i in 1:length(list_num_par)
        new_list = empty([Symbol(atomic_symbol(at[1])) => position(at[1])])
        pos = position(at, comp_split[i])
        symbols = [Symbol(atomic_symbol(at[i])) for i in comp_split[i]]
        for i in 1:length(pos)
            push!(new_list, symbols[i] => pos[i])
        end
        sys = atomic_system(new_list, cell_vectors(at), periodicity(at))
        components[i] = FastSystem(sys)
    end
    return components
end


"""
    split_components_by_chemical_species(at::AbstractSystem)

Split an `AbstractSystem` into multiple components based on the chemical species.

# Arguments
- `at::AbstractSystem`: The input `AbstractSystem` to be split.

# Returns
An array of `FastSystem` objects, each representing a component of the input system.

# Example
```jldoctest
julia> at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[2,1,3,4,5,6];particle_types=[:H,:O,:H,:Fe,:Au,:Cl])
FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF):
    bounding_box      : [ 5.94392        0        0;
                                0  5.94392        0;
                                0        0  5.94392]u"Å"

        .--------------.  
       /Au      Cl     |  
      / |HAu Fe        |  
     /  |     Cl Cl Cl |  
    *   |Cle           |  
    |   | Cl      H    |  
    |   |      OAuH    |  
    |FeFe-----------H--.  
    |  /          Au  /   
    | /      Au      /    
    |/              /     
    *--------------*      


julia> AbstractWalkers.split_components_by_chemical_species(at)
5-element Vector{FastSystem}:
 FastSystem(H₅, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å")
 FastSystem(O, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å")
 FastSystem(Cl₆, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å")
 FastSystem(Fe₄, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å")
 FastSystem(Au₅, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å")
```
"""
function split_components_by_chemical_species(at::AbstractSystem)
    list_species = [atomic_number(at, i) for i in 1:length(at)]
    species = sort!(unique(list_species))
    components = Array{FastSystem}(undef, length(species))
    for i in 1:length(species)
        new_list = empty([Symbol(atomic_symbol(at[1])) => position(at[1])])
        pos = position(at, findall(x->x==species[i],list_species))
        symbols = [Symbol(atomic_symbol(at[i])) for i in findall(x->x==species[i],list_species)]
        for i in 1:length(pos)
            push!(new_list, symbols[i] => pos[i])
        end
        sys = atomic_system(new_list, cell_vectors(at), periodicity(at))
        components[i] = FastSystem(sys)
    end
    return components
end

"""
    check_num_components(C::Int, list_num_par::Vector{Int}, frozen::Vector{Bool})

Check that the number of components matches the length of the list of number of particles and frozen particles.

# Arguments
- `C::Int`: The number of components.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

"""
function check_num_components(C::Int, list_num_par::Vector{Int}, frozen::Vector{Bool})
    if length(list_num_par) != C
        throw(ArgumentError("The number of components does not match the length of the list of number of particles."))
    elseif length(frozen) != C
        throw(ArgumentError("The number of components does not match the length of the list of frozen particles."))
    end
end

"""
    sort_components_by_atomic_number(at::AbstractSystem; merge_same_species=true)

Sorts the components of an `AbstractSystem` object `at` by their atomic number.

## Arguments
- `at::AbstractSystem`: The input `AbstractSystem` object.

## Keyword Arguments
- `merge_same_species::Bool=true`: Whether to merge components with the same species.

## Returns
- `list_num_par::Vector{Int64}`: A vector containing the number of each component species.
- `new_list::FastSystem`: A new `FastSystem` object with the sorted components.

The function first extracts the atomic numbers of the components in `at`. If `merge_same_species` is `true`, it sorts the unique species and counts the number of each species. If `merge_same_species` is `false`, it creates a list of species and their counts. It then sorts the species and counts by atomic number. Finally, it constructs a new `FastSystem` object with the sorted components and returns the list of species counts and the new `FastSystem` object.

# Examples
```jldoctest
julia> at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[2,1,3,4,5,6];particle_types=[:H,:O,:H,:Fe,:Au,:Cl])
FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF):
    bounding_box      : [ 5.94392        0        0;
                                0  5.94392        0;
                                0        0  5.94392]u"Å"

        .--------------.  
       /|     Cl    H  |  
      / |      Fe      |  
     /  |  Au     FeH  |  
    *   |   FeH  ACl   |  
    |   |    Cl     Au |  
    |   |            O |  
    |   .--Fe----------.  
    |  /H  Cl         /   
    | /          Au  /    
    |/Cl          Cl/     
    *--------------*      


julia> AbstractWalkers.sort_components_by_atomic_number(at; merge_same_species=false)
([2, 3, 1, 6, 4, 5], FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å"))

julia> AbstractWalkers.sort_components_by_atomic_number(at)
([5, 1, 6, 4, 5], FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å"))
```
"""
function sort_components_by_atomic_number(at::AbstractSystem; merge_same_species=true)
    list_species = [atomic_number(at, i) for i in 1:length(at)]
    new_list = empty([Symbol(atomic_symbol(at[1])) => position(at[1])])
    if merge_same_species
        species = sort!(unique(list_species))
        list_num_par = [count(x->x==s, list_species) for s in species]
    elseif !merge_same_species
        species = [list_species[1]; [list_species[i] for i in 2:length(list_species) if list_species[i] != list_species[i-1]]]
        list_num_par = Vector{Int64}()
        i = 1
        while i <= length(list_species)
            count = 1
            while i + count <= length(list_species) && list_species[i] == list_species[i + count]
                count += 1
            end
            push!(list_num_par, count)
            i += count
        end
        zipped = sort!(collect(zip(species, list_num_par)), by=first)
        list_num_par = [x[2] for x in zipped]
        species = unique([x[1] for x in zipped])
    end
    for i in 1:length(species)
        pos = position(at, findall(x->x==species[i],list_species))
        symbols = [Symbol(atomic_symbol(at[i])) for i in findall(x->x==species[i],list_species)]
        for i in 1:length(pos)
            push!(new_list, symbols[i] => pos[i])
        end
    end
    sys = atomic_system(new_list, cell_vectors(at), periodicity(at))
    return list_num_par, FastSystem(sys)
end

"""
    view_structure(at::AbstractSystem)
    view_structure(walker::AtomWalker)

Print an ASCII representation of the system.
"""
function view_structure(at::AbstractSystem)
    return println(visualize_ascii(at))
end

function view_structure(walker::AtomWalker)
    return println(visualize_ascii(walker.configuration))
end