# Atomistic walkers
abstract type AtomWalkers <: AbstractLiveSet end

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
- `LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::LennardJonesParametersSets; assign_energy=true)`: 
    Constructs a new `LJAtomWalkers` object with the given walkers and Lennard-Jones potential parameters. If `assign_energy=true`,
    the energy of each walker is assigned using the Lennard-Jones potential.

"""
struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    lj_potential::LennardJonesParametersSets
    function LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::LennardJonesParametersSets; assign_energy=true) where C
        if assign_energy
            [assign_energy!(walker, lj_potential) for walker in walkers]
        end
        return new(walkers, lj_potential)
    end
end

function Base.show(io::IO, walkers::LJAtomWalkers)
    println(io, "LJAtomWalkers($(eltype(walkers.walkers)), $(typeof(walkers.lj_potential))):")
    if length(walkers.walkers) > 10
        for i in 1:5
            println(io, "[$i] ", walkers.walkers[i])
        end
        println(io, "⋮\nOmitted ", length(walkers.walkers)-10, " walkers\n⋮\n")
        for i in length(walkers.walkers)-4:length(walkers.walkers)
            println(io, "[$i] ", walkers.walkers[i])
        end
    else
        for (ind, w) in enumerate(walkers.walkers)
            println(io, "[$ind] ", w)
        end
    end
    println(io, walkers.lj_potential)
end