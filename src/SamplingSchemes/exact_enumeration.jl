"""
    exact_enumeration(lattice::SLattice{G}, cutoff_radii::Tuple{Float64, Float64}, h::LatticeGasHamiltonian) where G

Enumerate all possible configurations of a lattice system and compute the energy of each configuration.

# Arguments
- `lattice::SLattice{G}`: The (starting) lattice system to enumerate. All possible configurations will be generated from this lattice system.
- `cutoff_radii::Tuple{Float64, Float64}`: The cutoff radii for the first and second nearest neighbors.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.

# Returns
- `energies::Vector{typeof(0.0u"eV")}`: An array of energies for each configuration.
- `configurations::Vector{SLattice{G}}`: An array of lattice system configurations for each configuration.
- `walkers::Vector{LatticeWalker}`: An array of lattice walkers for each configuration.
"""
function exact_enumeration(
    lattice::SLattice{G},
    h::ClassicalHamiltonian,
    ) where G

    primitive_lattice_vectors::Matrix{Float64} = lattice.lattice_vectors
    basis::Vector{Tuple{Float64, Float64, Float64}} = lattice.basis
    supercell_dimensions::Tuple{Int64, Int64, Int64} = lattice.supercell_dimensions
    periodicity::Tuple{Bool, Bool, Bool} = lattice.periodicity 
    adsorptions::Vector{Bool} = lattice.adsorptions
    cutoff_radii::Vector{Float64} = lattice.cutoff_radii
    number_occupied_sites::Int64 = sum(lattice.occupations)

    K, L, M = supercell_dimensions
    num_basis_sites = length(basis)
    total_sites = K * L * M * num_basis_sites

    # Generate all possible occupation configurations
    all_configs = collect(combinations(1:total_sites, number_occupied_sites))

    # Generate occupation vectors from configurations
    all_occupation_vectors = Vector{Vector{Bool}}()
    for config in all_configs
        occupations = falses(total_sites)
        occupations[config] .= true
        push!(all_occupation_vectors, Vector{Bool}(occupations))  # Convert BitVector to Vector{Bool}
    end

    # Generate SLattice objects for each configuration
    lattices = [SLattice{G}(primitive_lattice_vectors, basis, supercell_dimensions, periodicity, occupations, adsorptions, cutoff_radii) for occupations in all_occupation_vectors]

    # Generate LatticeWalker objects for each lattice system
    walkers = [LatticeWalker(lattice) for lattice in lattices]

    # Compute energies for each walker
    for walker in walkers
        e_interaction = interacting_energy(walker.configuration, h)
        walker.energy = e_interaction
    end

    # Extract energies and configurations
    energies = Array{typeof(0.0u"eV")}(undef, length(walkers))
    configurations = Array{SLattice{G}}(undef, length(walkers))

    for (i, walker) in enumerate(walkers)
        energies[i] = walker.energy
        configurations[i] = walker.configuration
    end

    return energies, configurations, walkers
end