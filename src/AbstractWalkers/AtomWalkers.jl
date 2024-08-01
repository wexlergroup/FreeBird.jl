"""
    mutable struct AtomWalker

The `AtomWalker` struct represents a walker composed of atoms/molecules.

# Fields
- `configuration::FastSystem`: The configuration of the walker.
- `energy::typeof(0.0u"eV")`: The energy of the walker.
- `iter::Int64`: The current iteration number of the walker.
- `list_num_par::Vector{Int64}`: The list of the number of particles for each component.
- `frozen::Vector{Bool}`: A boolean vector indicating whether each component is frozen or not.
- `energy_frozen_part::typeof(0.0u"eV")`: The energy of the frozen particles in the walker, serves as a constant energy offset
    to the interacting part of the system.

# Constructor
- `AtomWalker(configuration::FastSystem; energy=0.0u"eV", iter=0, list_num_par=zeros(Int,C), frozen=zeros(Bool,C), energy_frozen_part=0.0u"eV")`: 
    Constructs a new `AtomWalker` object with the given configuration and optional parameters.

"""
mutable struct AtomWalker{C} <: AbstractWalker
    configuration::FastSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    list_num_par::Vector{Int64}
    frozen::Vector{Bool}
    energy_frozen_part::typeof(0.0u"eV")
    function AtomWalker{C}(configuration::FastSystem; energy=0.0u"eV", iter=0, list_num_par=zeros(Int,C), frozen=zeros(Bool,C), energy_frozen_part=0.0u"eV") where C
        if length(list_num_par) != C
            throw(ArgumentError("the length of the list_num_par is not compatible with the number of components."))
        elseif length(frozen) != C
            throw(ArgumentError("the length of the frozen is not compatible with the number of components."))
        end
        if C==1
            list_num_par = [length(configuration)]
        end
        if sum(list_num_par) != length(configuration)
            throw(ArgumentError("the sum of the list_num_par is not compatible with the length of the configuration."))
        end
        return new{C}(configuration, energy, iter, list_num_par, frozen, energy_frozen_part)
    end
end

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
    AtomWalker(configuration::FastSystem; freeze_species::Vector{Symbol}=Symbol[], merge_same_species=true)

Constructs an `AtomWalker` object with the given configuration.

# Arguments
- `configuration::FastSystem`: The configuration of the walker.
- `freeze_species::Vector{Symbol}`: A vector of species to freeze.
- `merge_same_species::Bool`: A boolean indicating whether to merge the same species into one component.

# Returns
- `AtomWalker{C}`: The constructed `AtomWalker` object.

# Example
```jldoctest
julia> at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[2,1,3,4,5,6];particle_types=[:H,:O,:H,:Fe,:Au,:Cl])
FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF):
    bounding_box      : [ 5.94392        0        0;
                                0  5.94392        0;
                                0        0  5.94392]u"Å"

        .--------------.  
       /|Fel           |  
      / H   H   Cl     |  
     /  Hu   O         |  
    *   |       Au   Fe|  
    |   |FeCl        Fe|  
    |   |        Au    |  
    |   .---------Au---.  
    |  /           H  /   
    | Au Cl          /    
    |/              /     
    *--------------*      

julia> AtomWalker(at;freeze_species=[:H],merge_same_species=false)
AtomWalker{6}(FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å"), 0.0 eV, 0, [2, 3, 1, 6, 4, 5], Bool[1, 1, 0, 0, 0, 0], 0.0 eV)

julia> AtomWalker(at;freeze_species=[:H],merge_same_species=true)
AtomWalker{5}(FastSystem(Au₅Cl₆Fe₄H₅O, periodic = FFF, bounding_box = [[5.943921952763129, 0.0, 0.0], [0.0, 5.943921952763129, 0.0], [0.0, 0.0, 5.943921952763129]]u"Å"), 0.0 eV, 0, [5, 1, 6, 4, 5], Bool[1, 0, 0, 0, 0], 0.0 eV)
```

"""
function AtomWalker(configuration::FastSystem; freeze_species::Vector{Symbol}=Symbol[], merge_same_species=true)
    list_num_par, configuration = sort_components_by_atomic_number(configuration, merge_same_species=merge_same_species)
    C = length(list_num_par)
    frozen = zeros(Bool, C)
    if !isempty(freeze_species)
        elements = unique(atomic_symbol(configuration))
        @show elements
        for species in freeze_species
            if !(species in elements)
                throw(ArgumentError("The species $species is not in the system."))
            end
        end
        components = split_components(configuration, list_num_par)
        for i in 1:length(freeze_species)
            for j in 1:C
                if freeze_species[i] in atomic_symbol(components[j])
                    frozen[j] = true
                end
            end
        end
    end
    return AtomWalker{C}(configuration, list_num_par=list_num_par, frozen=frozen)
end



"""
    update_walker!(walker::AtomWalker, key::Symbol, value)

Update the properties of an AtomWalker object.

A convenient function that updates the value of a specific property of an AtomWalker object.

# Arguments
- `walker::AtomWalker`: The AtomWalker object to be updated.
- `key::Symbol`: The key of the property to be updated.
- `value`: The new value of the property.

# Returns
- `walker::AtomWalker`: The updated AtomWalker object.

# Example
```julia
update_walker!(walker, :energy, 10.0u"eV")
update_walker!(walker, :iter, 1)
```

"""
function update_walker!(walker::AtomWalker, key::Symbol, value)
    setproperty!(walker, key, value)
    return walker
end