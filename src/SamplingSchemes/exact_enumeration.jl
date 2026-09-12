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
    sort!(comp_list)
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
    exact_enumeration(lattice::SLattice{G}, cutoff_radii::Tuple{Float64, Float64}, h::ClassicalHamiltonian) where G

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


"""
    grand_canonical_exact_enumeration(lattice::MLattice{1,G}, h::ClassicalHamiltonian;
                                      max_sites::Int=20) where G

Enumerate **all** \$2^M\$ occupations of a single-component lattice and return
`(energies, numbers)`: the bare energy and the particle number of every
microstate, in mask order.

This is the grand-canonical counterpart to [`exact_enumeration`](@ref), and
exists because that function cannot serve as a grand-canonical reference:
`enumerate_lattices(::SLattice)` is **fixed-N** — it takes
`number_occupied_sites = sum(lattice.components[1])` and enumerates
`combinations(1:total_sites, number_occupied_sites)` — so it covers one particle
number, not the full state space a grand-canonical ensemble samples over.

Nothing about \$\\mu\$ or \$\\beta\$ is baked in. Both are properties of the
ensemble you evaluate, not of the spectrum, so the caller supplies them:

```julia
E, N = grand_canonical_exact_enumeration(lattice, h)
Ω = ustrip.(E) .- μ .* N
w = exp.(-β .* Ω)
mean_N = sum(w .* N) / sum(w)
```

That separation is the point of promoting this out of the test files: the same
spectrum serves every \$(\\mu, \\beta)\$ a test or a regression wants to check,
and there is one implementation of it to be wrong rather than four.

# Arguments
- `lattice::MLattice{1,G}`: The lattice whose site count and geometry define the
  state space. Its current occupation is irrelevant — every occupation is
  visited — and it is not mutated.
- `h::ClassicalHamiltonian`: The Hamiltonian used for every microstate.
- `max_sites::Int=20`: Guard against accidentally asking for \$2^{40}\$ states.
  Raise it deliberately if you mean it.

# Returns
- `energies::Vector`: Bare \$E\$ (Unitful, as `interacting_energy` returns it)
  for each of the \$2^M\$ microstates.
- `numbers::Vector{Int}`: \$N\$ for the same microstates, in the same order.

Mask order: microstate `i` (1-based) has site `j` occupied iff bit `j-1` of
`i-1` is set, so `numbers[i] == count_ones(i-1)`.
"""
function grand_canonical_exact_enumeration(lattice::MLattice{1,G},
                                           h::ClassicalHamiltonian;
                                           max_sites::Int=20) where G
    M = num_sites(lattice)
    if M > max_sites
        throw(ArgumentError(
            "grand_canonical_exact_enumeration: $M sites means 2^$M microstates. " *
            "Raise max_sites (currently $max_sites) if that is really intended."))
    end

    n_states = 1 << M
    template = deepcopy(lattice)
    energies = Vector{typeof(interacting_energy(template, h))}(undef, n_states)
    numbers = Vector{Int}(undef, n_states)

    # One lattice per task, strided over the masks. A task-local copy rather
    # than an array indexed by threadid(), which is not safe under task
    # migration; and mutating one lattice per task rather than materialising
    # 2^M of them, which is what the hand-rolled loops in the test files did.
    nchunks = max(1, Threads.nthreads())
    Threads.@threads for c in 1:nchunks
        lat = deepcopy(lattice)
        occ = lat.components[1]
        for mask in (c - 1):nchunks:(n_states - 1)
            for site in 1:M
                occ[site] = ((mask >> (site - 1)) & 1) == 1
            end
            energies[mask + 1] = interacting_energy(lat, h)
            numbers[mask + 1] = count_ones(mask)
        end
    end

    return energies, numbers
end

function grand_canonical_exact_enumeration(lattice::MLattice{C,G},
                                           ::ClassicalHamiltonian;
                                           max_sites::Int=20) where {C,G}
    throw(ArgumentError(
        "grand_canonical_exact_enumeration is defined for single-component " *
        "lattices; got MLattice{$C,$G}. The grand-canonical samplers are " *
        "single-component too (they read components[1])."))
end
