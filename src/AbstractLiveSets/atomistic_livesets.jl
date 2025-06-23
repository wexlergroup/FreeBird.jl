# Atomistic walkers
abstract type AtomWalkers <: AbstractLiveSet end

"""
    assign_energy!(walker::AtomWalker, lj::PotentialParameterSets)

Assigns the energy to the given `walker` using the Lennard-Jones parameters `lj`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `lj::PotentialParameterSets`: The Lennard-Jones parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_energy!(walker::AtomWalker, lj::PotentialParameterSets)
    # walker.energy_frozen_part = frozen_energy(walker.configuration, lj, walker.list_num_par, walker.frozen)
    if lj isa SMD_LJParameters
        walker.energy = interacting_energy(walker.configuration, lj)
    else
        walker.energy = interacting_energy(walker.configuration, lj, walker.list_num_par, walker.frozen) + walker.energy_frozen_part
    end
    return walker
end

"""
    assign_frozen_energy!(walker::AtomWalker, lj::PotentialParameterSets)

Assigns the frozen energy to the given `walker` using the Lennard-Jones parameters `lj`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `lj::PotentialParameterSets`: The Lennard-Jones parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_frozen_energy!(walker::AtomWalker, lj::PotentialParameterSets)
    walker.energy_frozen_part = frozen_energy(walker.configuration, lj, walker.list_num_par, walker.frozen)
    return walker
end


"""
    struct LJAtomWalkers <: AtomWalkers

The `LJAtomWalkers` struct represents a collection of atom walkers that interact with each other using the Lennard-Jones potential.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `lj_potential::PotentialParameterSets`: The Lennard-Jones potential parameters. See `PotentialParameterSets`.

# Constructor
- `LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::PotentialParameterSets; assign_energy=true)`: 
    Constructs a new `LJAtomWalkers` object with the given walkers and Lennard-Jones potential parameters. If `assign_energy=true`,
    the energy of each walker is assigned using the Lennard-Jones potential.

"""
struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    lj_potential::PotentialParameterSets
    function LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::PotentialParameterSets; assign_energy=true, const_frozen_part=true) where C
        if lj_potential isa SMD_LJParameters
            frozen_part_energy = 0.0 * unit(lj_potential.epsilon)
        elseif const_frozen_part && !isempty(walkers) 
            frozen_part_energy = frozen_energy(walkers[1].configuration, lj_potential, walkers[1].list_num_par, walkers[1].frozen)
        end
        if assign_energy
            for walker in walkers
                if const_frozen_part
                    walker.energy_frozen_part = frozen_part_energy
                else
                    assign_frozen_energy!(walker, lj_potential)
                end
                assign_energy!(walker, lj_potential) # comes after assign_frozen_energy!
            end
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