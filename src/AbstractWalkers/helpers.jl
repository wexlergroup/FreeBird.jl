# convenience functions for handling atomistic walkers

"""
    split_components(at::AbstractSystem, list_num_par::Vector{Int})

Split the system into components based on the number of particles in each component.

# Arguments
- `at::AbstractSystem`: The system to split.
- `list_num_par::Vector{Int}`: The number of particles in each component.

# Returns
- `components`: An array of `FastSystem` objects representing the components of the system.

An empty system, or a zero-count component, yields a zero-atom `FastSystem` carrying the
parent system's cell and periodicity.
"""
function split_components(at::AbstractSystem, list_num_par::Vector{Int})
    components = Array{FastSystem}(undef, length(list_num_par))
    comp_cut = vcat([0],cumsum(list_num_par))
    comp_split = [comp_cut[i]+1:comp_cut[i+1] for i in 1:length(list_num_par)]
    for i in 1:length(list_num_par)
        new_list = Pair{Symbol, eltype(position(at, :))}[]
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
An empty system returns an empty array.

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
        new_list = Pair{Symbol, eltype(position(at, :))}[]
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

The function first extracts the atomic numbers of the components in `at`. If `merge_same_species` is `true`, it sorts the unique species and counts the number of each species. If `merge_same_species` is `false`, it creates a list of species and their counts. It then sorts the species and counts by atomic number. Finally, it constructs a new `FastSystem` object with the sorted components and returns the list of species counts and the new `FastSystem` object. An empty system returns an empty `list_num_par` and a zero-atom system, under either flag.

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
    new_list = Pair{Symbol, eltype(position(at, :))}[]
    if merge_same_species
        species = sort!(unique(list_species))
        list_num_par = [count(x->x==s, list_species) for s in species]
    elseif !merge_same_species
        species = isempty(list_species) ? empty(list_species) : [list_species[1]; [list_species[i] for i in 2:length(list_species) if list_species[i] != list_species[i-1]]]
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
    insert_particle!(walker::AtomWalker{1}, pos, species)

Append one particle of `species` (a `Symbol` or `ChemicalSpecies`) at position `pos` to the
walker's configuration, updating `list_num_par` in the same call. The two mutations are kept
in lockstep here so the configuration arrays and the particle count cannot drift apart. The
operation is purely structural: the walker's `energy` is not touched and the component's
frozen flag is not consulted; callers own the energy bookkeeping (see
`MC_grand_canonical_walk!`).

# Arguments
- `walker::AtomWalker{1}`: The single-component walker to extend.
- `pos`: The new particle's position (any 3-vector of lengths accepted by the configuration).
- `species`: The chemical identity of the new particle.

# Returns
- `walker::AtomWalker{1}`: The updated walker.
"""
function insert_particle!(walker::AtomWalker{1}, pos, species)
    config = walker.configuration
    sp = species isa ChemicalSpecies ? species : ChemicalSpecies(species)
    push!(config.position, pos)
    push!(config.species, sp)
    push!(config.mass, mass(sp))
    walker.list_num_par[1] += 1
    return walker
end

"""
    remove_particle!(walker::AtomWalker{1}, i::Int)

Remove particle `i` from the walker's configuration, order-preserving, updating
`list_num_par` in the same call. Purely structural, like `insert_particle!`: the walker's
`energy` is not touched.

# Arguments
- `walker::AtomWalker{1}`: The single-component walker to shrink.
- `i::Int`: The index of the particle to remove.

# Returns
- `walker::AtomWalker{1}`: The updated walker.
"""
function remove_particle!(walker::AtomWalker{1}, i::Int)
    config = walker.configuration
    checkbounds(config.position, i)
    deleteat!(config.position, i)
    deleteat!(config.species, i)
    deleteat!(config.mass, i)
    walker.list_num_par[1] -= 1
    return walker
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