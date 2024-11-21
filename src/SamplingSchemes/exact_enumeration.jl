"""
    exact_enumeration(lattice::SLattice{G}, cutoff_radii::Tuple{Float64, Float64}, h::LatticeGasHamiltonian) where G

Enumerate all possible configurations of a lattice system and compute the energy of each configuration.

# Arguments
- `lattice::SLattice{G}`: The (starting) lattice system to enumerate. All possible configurations will be generated from this lattice system.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.

# Returns
- `DataFrame`: A DataFrame containing the energy and configuration of each configuration.
- `LatticeGasWalkers`: A collection of lattice walkers for each configuration.
"""
function exact_enumeration(lattice::SLattice{G}, h::ClassicalHamiltonian) where G

    number_occupied_sites::Int64 = sum(lattice.occupations)
    total_sites = length(lattice.basis) * prod(lattice.supercell_dimensions)

    # Generate all possible occupation configurations
    all_configs = combinations(1:total_sites, number_occupied_sites)

    # Generate occupation vectors from configurations
    lattices = [deepcopy(lattice) for _ in 1:length(all_configs)]

    for (ind, config) in enumerate(all_configs)
        occupations = lattices[ind].occupations
        occupations[config] .= true
    end

    ls = LatticeGasWalkers(LatticeWalker.(lattices), h)

    # Extract energies and configurations
    energies = [wk.energy for wk in ls.walkers]
    configurations = [wk.configuration.occupations for wk in ls.walkers]

    df = DataFrame()
    df.energy = energies
    df.config = configurations

    return df, ls
end