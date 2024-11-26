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

    # flush occupancy
    lattice.occupations .= false

    # Generate occupation vectors from configurations
    lattices = [deepcopy(lattice) for _ in 1:length(all_configs)]

    Threads.@threads for (ind, config) in collect(enumerate(all_configs))
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

# recursive function to generate all unique permutations of a list
# https://stackoverflow.com/questions/65051953/julia-generate-all-non-repeating-permutations-in-set-with-duplicates
function unique_permutations(x::T, prefix=T()) where T
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i-1]
                continue
            end
            append!(t, unique_permutations([x[begin:i-1];x[i+1:end]], [prefix; x[i]]))
        end
        return t
    end
end

function exact_enumeration(lattice::MLattice{C,G}, h::ClassicalHamiltonian) where {C,G}

    total_sites = length(lattice.basis) * prod(lattice.supercell_dimensions)

    # setup a vector of all components
    comp_list = zeros(Int, total_sites)
    for i in 1:C
        comp_list += lattice.components[i] * i
    end
    # Generate all possible occupation configurations
    all_configs = unique_permutations(comp_list)

    # flush occupancy
    for i in 1:C
        lattice.components[i].= false
    end

    # Generate occupation vectors from configurations
    lattices = [deepcopy(lattice) for _ in 1:length(all_configs)]

    Threads.@threads for (ind, config) in collect(enumerate(all_configs))
        for i in 1:C
            lattices[ind].components[i] = config .== i
        end
    end

    ls = LatticeGasWalkers(LatticeWalker.(lattices), h)

    # Extract energies and configurations
    energies = [wk.energy for wk in ls.walkers]
    configurations = [wk.configuration.components for wk in ls.walkers]

    df = DataFrame()
    df.energy = energies
    df.config = configurations

    return df, ls
end