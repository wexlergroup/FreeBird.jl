# Lattice walkers
abstract type LatticeWalkers <: AbstractLiveSet end


"""
    assign_energy!(walker::LatticeWalker{C}, hamiltonian::LatticeGasHamiltonian; perturb_energy::Float64=0.0)

Assigns energy to the given `walker` based on the `hamiltonian`. If `perturb_energy` is non-zero, a small random perturbation is added to the energy.

# Arguments
- `walker::LatticeWalker{C}`: The walker to assign energy to.
- `hamiltonian::LatticeGasHamiltonian`: The Hamiltonian used to calculate the energy.
- `perturb_energy::Float64=0.0`: The amount of random perturbation to add to the energy.

# Returns
- `walker::LatticeWalker{C}`: The walker with the assigned energy.
"""
function assign_energy!(walker::LatticeWalker{C}, hamiltonian::ClassicalHamiltonian; perturb_energy::Float64=0.0) where C
    # Assign the energy to the walker and, if perturb_energy is non-zero, give all walkers a small random (positive or negative) perturbation
    walker.energy = interacting_energy(walker.configuration, hamiltonian) + perturb_energy * (rand() - 0.5) * unit(walker.energy)
    return walker
end

"""
    struct LatticeGasWalkers <: LatticeWalkers

The `LatticeGasWalkers` struct represents a collection of lattice walkers for a lattice gas system. It is a subtype of `LatticeWalkers`.

# Fields
- `walkers::Vector{LatticeWalker{C}}`: A vector of lattice walkers.
- `hamiltonian::LatticeGasHamiltonian`: The lattice gas Hamiltonian associated with the walkers.

# Constructors
- `LatticeGasWalkers(walkers::Vector{LatticeWalker{C}}, hamiltonian::LatticeGasHamiltonian; assign_energy=true, perturb_energy::Float64=0.0)`: Constructs a new `LatticeGasWalkers` object with the given walkers and Hamiltonian. If `assign_energy` is `true`, the energy of each walker is assigned using the provided Hamiltonian. The optional `perturb_energy` parameter can be used to add a small perturbation to the assigned energy.

"""
struct  LatticeGasWalkers <: LatticeWalkers
    walkers::Vector{LatticeWalker{C}} where C
    hamiltonian::ClassicalHamiltonian
    function LatticeGasWalkers(walkers::Vector{LatticeWalker{C}}, hamiltonian::ClassicalHamiltonian; assign_energy=true, perturb_energy::Float64=0.0) where C
        if assign_energy
            [assign_energy!(walker, hamiltonian; perturb_energy=perturb_energy) for walker in walkers]
        end
        return new(walkers, hamiltonian)
    end
end

function Base.show(io::IO, ls::LatticeGasWalkers)
    println(io, "LatticeGasWalkers($(eltype(ls.walkers)), $(typeof(ls.hamiltonian))):")
    println(io, "    lattice_vectors:      ", ls.walkers[1].configuration.lattice_vectors)
    println(io, "    supercell_dimensions: ", ls.walkers[1].configuration.supercell_dimensions)
    println(io, "    periodicity:          ", ls.walkers[1].configuration.periodicity)
    println(io, "    basis:                ", ls.walkers[1].configuration.basis)
    if length(ls.walkers) > 10
        for i in 1:5
            println(io, "[$i] ", "energy = ", ls.walkers[i].energy, ", iter = ", ls.walkers[i].iter)
            print_lattice_walker_in_walkers(io, ls.walkers[i])
        end
        println(io, "⋮\nOmitted ", length(ls.walkers)-10, " walkers\n⋮\n")
        for i in length(ls.walkers)-4:length(ls.walkers)
            println(io, "[$i] ", "energy = ", ls.walkers[i].energy, ", iter = ", ls.walkers[i].iter)
            print_lattice_walker_in_walkers(io, ls.walkers[i])
        end
    else
        for (i, walker) in enumerate(ls.walkers)
            println(io, "[$i] ", "energy = ", ls.walkers[i].energy, ", iter = ", ls.walkers[i].iter)
            print_lattice_walker_in_walkers(io, ls.walkers[i])
        end
    end
    println(io)
    println(io, ls.hamiltonian)
end

function print_lattice_walker_in_walkers(io::IO, walker::LatticeWalker{C}) where C
    println(io, "    occupations:")
    if walker.configuration isa MLattice
        # AbstractWalkers.print_lattice(io, walker.configuration, walker.configuration.components)
        println(io, "      components:")
        for (i, component) in enumerate(walker.configuration.components)
            println(io, "        component $i:")
            println(io, "          ", component)
        end
    else
        AbstractWalkers.print_lattice(io, walker.configuration, walker.configuration.occupations)
    end
end

function Base.show(io::IO, list::Vector{LatticeWalker{C}}) where C
    println(io, "Vector{LatticeWalker{$C}}(", length(list), "):")
    if length(list) > 10
        for i in 1:5
            println(io, "[$i] ", "energy = ", list[i].energy, ", iter = ", list[i].iter)
            AbstractWalkers.print_occupation(io, list[i].configuration)
        end
        println(io, "⋮\nOmitted ", length(list)-10, " walkers\n⋮\n")
        for i in length(list)-4:length(list)
            println(io, "[$i] ", "energy = ", list[i].energy, ", iter = ", list[i].iter)
            AbstractWalkers.print_occupation(io, list[i].configuration)
        end
    else
        for (i, walker) in enumerate(list)
            println(io, "[$i] ", "energy = ", list[i].energy, ", iter = ", list[i].iter)
            AbstractWalkers.print_occupation(io, list[i].configuration)
        end
    end
end