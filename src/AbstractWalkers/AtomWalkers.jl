abstract type AtomWalkers end

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
mutable struct AtomWalker{C}
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


"""
    assign_energy!(walker::AtomWalker, lj::LennardJonesParametersSets)

Assigns the energy to the given `walker` using the Lennard-Jones parameters `lj`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `lj::LennardJonesParametersSets`: The Lennard-Jones parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_energy!(walker::AtomWalker, lj::LennardJonesParametersSets)
    walker.energy_frozen_part = frozen_energy(walker.configuration, lj, walker.list_num_par, walker.frozen)
    walker.energy = interacting_energy(walker.configuration, lj, walker.list_num_par, walker.frozen) + walker.energy_frozen_part
    return walker
end


"""
    struct LJAtomWalkers <: AtomWalkers

The `LJAtomWalkers` struct represents a collection of atom walkers that interact with each other using the Lennard-Jones potential.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `lj_potential::LennardJonesParametersSets`: The Lennard-Jones potential parameters. See `LennardJonesParametersSets`.

# Constructor
- `LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::LennardJonesParametersSets) where C`: Constructs a new `LJAtomWalkers` object with the given walkers and Lennard-Jones potential parameters.

"""
struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    lj_potential::LennardJonesParametersSets
    function LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::LennardJonesParametersSets) where C
        [assign_energy!(walker, lj_potential) for walker in walkers]
        return new(walkers, lj_potential)
    end
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