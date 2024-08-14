# Lattice walkers
abstract type LatticeWalkers <: AbstractLiveSet end

function assign_energy!(walker::LatticeWalker, hamiltonian::LatticeGasHamiltonian; perturb_energy::Float64=0.0)
    walker.energy = interacting_energy(walker.configuration, hamiltonian) * (1 + perturb_energy * (round(rand(), sigdigits=1) - 0.5))
    return walker
end

struct  LatticeGasWalkers <: LatticeWalkers
    walkers::Vector{LatticeWalker}
    hamiltonian::LatticeGasHamiltonian
    function LatticeGasWalkers(walkers::Vector{LatticeWalker}, hamiltonian::LatticeGasHamiltonian; assign_energy=true, perturb_energy::Float64=0.0)
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

function print_lattice_walker_in_walkers(io::IO, walkers::LatticeWalker)
    println(io, "    occupations:")
    AbstractWalkers.print_lattice(io, walkers.configuration, walkers.configuration.occupations)
    # println(io, "    adsorptions:")
    # AbstractWalkers.print_lattice(io, walkers.configuration, walkers.configuration.adsorptions)
end