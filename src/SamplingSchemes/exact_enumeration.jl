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

function enumerate_lattices(init_lattice::MLattice{C,G}) where {C,G}
    total_sites = length(init_lattice.basis) * prod(init_lattice.supercell_dimensions)

    # setup a vector of all components
    comp_list = zeros(Int, total_sites)
    for i in 1:C
        comp_list += init_lattice.components[i] * i
    end
    # Generate all possible occupation configurations
    all_configs = unique_permutations(comp_list)

    lattice = deepcopy(init_lattice)

    # flush occupancy
    for i in 1:C
        lattice.components[i].= false
    end

    # Generate occupation vectors from configurations
    lattices = [deepcopy(lattice) for _ in eachindex(all_configs)]

    Threads.@threads for (ind, config) in collect(enumerate(all_configs))
        for i in 1:C
            lattices[ind].components[i] = config .== i
        end
    end

    return lattices
end


function enumerate_lattices(init_lattice::SLattice{G}) where {G}
    
    # @debug "SLattice routine called"

    lattice = deepcopy(init_lattice)

    number_occupied_sites::Int64 = sum(lattice.components[1])
    total_sites = length(lattice.basis) * prod(lattice.supercell_dimensions)

    # Generate all possible occupation configurations
    # 3X faster than using unique_permutations for SLattice
    all_configs = combinations(1:total_sites, number_occupied_sites)

    # flush occupancy
    # lattice.components[1] .= false

    # Generate occupation vectors from configurations
    lattices = Vector{typeof(lattice)}(undef, length(all_configs))

    # non-threaded version - debug purposes only
    # for (ind, config) in collect(enumerate(all_configs))
    #     lattice.components[1] .= false # flush occupancy
    #     lattice.components[1][config] .= true
    #     lattices[ind] = deepcopy(lattice)
    # end

    all_configs = collect(all_configs)
    
    Threads.@threads for ind in eachindex(all_configs)
        new_lattice = deepcopy(lattice)
        new_lattice.components[1] .= false # flush occupancy
        new_lattice.components[1][all_configs[ind]] .= true
        # lattice.components[1][config] .= true
        lattices[ind] = new_lattice
    end

    return lattices
end

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
function exact_enumeration(lattice::MLattice{C,G}, h::ClassicalHamiltonian) where {C,G}

    lattices = enumerate_lattices(lattice)

    ls = LatticeGasWalkers(LatticeWalker.(lattices), h)

    # Extract energies and configurations
    energies = Vector{typeof(ls.walkers[1].energy)}(undef, length(ls.walkers))
    configurations = Vector{Vector{Vector{Bool}}}(undef, length(ls.walkers))

    Threads.@threads for i in eachindex(ls.walkers)
        energies[i] = ls.walkers[i].energy
        configurations[i] = ls.walkers[i].configuration.components
    end

    df = DataFrame()
    df.energy = energies
    df.config = configurations

    return df, ls
end


