abstract type LatticeWalkers end

"""
compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, positions::Matrix{Float64}, cutoff_radii::Tuple{Float64, Float64}, periodicity::Vector{Bool})

Compute the nearest and next-nearest neighbors for each atom in a 3D lattice.

# Arguments
- `supercell_lattice_vectors::Matrix{Float64}`: The lattice vectors of the supercell.
- `positions::Matrix{Float64}`: The positions of the atoms in the supercell.
- `periodicity::Vector{Bool}`: A Boolean vector of length three indicating periodicity in each dimension (true for periodic, false for non-periodic).
- `cutoff_radii::Tuple{Float64, Float64}`: The cutoff radii for the first and second nearest neighbors.

# Returns
- `neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}`: A vector of tuples containing the indices of the first and second nearest neighbors for each atom.

"""

function compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, positions::Matrix{Float64}, periodicity::Tuple{Bool, Bool, Bool}, cutoff_radii::Tuple{Float64, Float64})
    neighbors = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, size(positions, 1))
    num_atoms = size(positions, 1)
    
    # Compute reciprocal lattice vectors for minimum image convention
    a1 = supercell_lattice_vectors[:, 1]
    a2 = supercell_lattice_vectors[:, 2]
    a3 = supercell_lattice_vectors[:, 3]
    reciprocal_lattice_vectors = inv([a1 a2 a3])

    first_nearest_distance = cutoff_radii[1]
    second_nearest_distance = cutoff_radii[2]

    for i in 1:num_atoms
        first_neighbors = Int[]
        second_neighbors = Int[]
        pos_i = positions[i, :]
        
        for j in 1:num_atoms
            if i != j
                pos_j = positions[j, :]
                dx = pos_j[1] - pos_i[1]
                dy = pos_j[2] - pos_i[2]
                dz = pos_j[3] - pos_i[3]

                # Apply minimum image convention using reciprocal lattice vectors
                dr = [dx, dy, dz]
                fractional_dr = reciprocal_lattice_vectors * dr
                
                for k in 1:3
                    if periodicity[k]
                        fractional_dr[k] -= round(fractional_dr[k])
                    end
                end
                
                dr = supercell_lattice_vectors * fractional_dr

                distance = norm(dr)
                
                if distance <= first_nearest_distance
                    push!(first_neighbors, j)
                elseif distance <= second_nearest_distance
                    push!(second_neighbors, j)
                end
            end
        end
        
        neighbors[i] = (first_neighbors, second_neighbors)
    end
    
    return neighbors
end

"""
mutable struct LatticeSystem

The `LatticeSystem` struct represents a 3D lattice system.

# Fields
- `lattice_vectors::Matrix{Float64}`: The lattice vectors of the system.
- `positions::Matrix{Float64}`: The positions of the atoms in the system.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `occupations::Vector{Bool}`: A vector of booleans indicating whether each site is occupied.
- `neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}`: A vector of tuples containing the indices of the first and second nearest neighbors for each atom.
- `adsorptions::Vector{Bool}`: A vector of booleans indicating whether each site is adsorbed.

# Constructors
```julia
LatticeSystem(lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64, Int64}, occupations::Vector{Bool}, adsorptions::Vector{Bool}, cutoff_radii::Tuple{Float64, Float64}, periodicity::Vector{Bool})
```
Create a new `LatticeSystem` with the given lattice vectors, basis, supercell dimensions, occupations, adsorptions, cutoff radii, and periodicity.

"""

mutable struct LatticeSystem
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    occupations::Vector{Bool}
    neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}
    adsorptions::Vector{Bool}

    function LatticeSystem(
        lattice_vectors::Matrix{Float64},
        basis::Vector{Tuple{Float64, Float64, Float64}},
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        periodicity::Tuple{Bool, Bool, Bool},
        occupations::Vector{Bool},
        adsorptions::Vector{Bool},
        cutoff_radii::Tuple{Float64, Float64}
    )
        num_basis_sites = length(basis)
        num_supercell_sites = supercell_dimensions[1] * supercell_dimensions[2] * supercell_dimensions[3] * num_basis_sites
        
        positions = zeros(Float64, num_supercell_sites, 3)
        index = 1

        a1 = lattice_vectors[:, 1]
        a2 = lattice_vectors[:, 2]
        a3 = lattice_vectors[:, 3]

        for k in 1:supercell_dimensions[3]
            for j in 1:supercell_dimensions[2]
                for i in 1:supercell_dimensions[1]
                    for (bx, by, bz) in basis
                        x = (i - 1) * a1[1] + (j - 1) * a2[1] + (k - 1) * a3[1] + bx
                        y = (i - 1) * a1[2] + (j - 1) * a2[2] + (k - 1) * a3[2] + by
                        z = (i - 1) * a1[3] + (j - 1) * a2[3] + (k - 1) * a3[3] + bz
                        positions[index, :] = [x, y, z]
                        index += 1
                    end
                end
            end
        end

        if length(occupations) != size(positions, 1)
            throw(ArgumentError("Length of occupations vector must match the number of lattice sites"))
        end

        if length(adsorptions) != size(positions, 1)
            throw(ArgumentError("Length of adsorptions vector must match the number of lattice sites"))
        end

        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, positions, periodicity, cutoff_radii)
        
        return new(lattice_vectors, positions, supercell_dimensions, occupations, neighbors, adsorptions)
    end
end

"""
    mutable struct LatticeWalker

The `LatticeWalker` struct represents a walker on a 3D lattice.

# Fields
- `configuration::LatticeSystem`: The configuration of the walker.
- `energy::Float64`: The energy of the walker.
- `iter::Int64`: The current iteration number of the walker.

# Constructor
```julia
LatticeWalker(configuration::LatticeSystem; energy=0.0, iter=0)
```
Create a new `LatticeWalker` with the given configuration and optional energy and iteration number.

"""

mutable struct LatticeWalker
    configuration::LatticeSystem
    energy::Float64
    iter::Int64
    function LatticeWalker(configuration::LatticeSystem; energy=0.0, iter=0)
        return new(configuration, energy, iter)
    end
end

"""
    interaction_energy(at::LatticeSystem, lg::LGHamiltonian)

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `at::LatticeSystem`: The lattice configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""

function interaction_energy(at::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64)
    e_adsorption = sum(at.occupations .& at.adsorptions) * adsorption_energy
    e_nn = 0.0
    e_nnn = 0.0

    for index in 1:length(at.occupations)
        if at.occupations[index]
            # Compute nearest-neighbor interaction energy
            for nn in at.neighbors[index][1]
                if at.occupations[nn]
                    e_nn += nn_energy / 2
                end
            end

            # Compute next-nearest-neighbor interaction energy
            for nnn in at.neighbors[index][2]
                if at.occupations[nnn]
                    e_nnn += nnn_energy / 2
                end
            end
        end
    end

    e_interaction = e_adsorption + e_nn + e_nnn
    return e_interaction
end

"""
    exact_enumeration(primitive_lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64}, number_occupied_sites::Int64, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, cutoff_radii::Tuple{Float64, Float64})

Enumerate all possible configurations of a lattice system and compute the energy of each configuration.

# Arguments
- `primitive_lattice_vectors::Matrix{Float64}`: The primitive lattice vectors of the system.
- `basis::Vector{Tuple{Float64, Float64}}`: The basis of the system.
- `supercell_dimensions::Tuple{Int64, Int64}`: The dimensions of the supercell.
- `number_occupied_sites::Int64`: The number of occupied sites in each configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.
- `cutoff_radii::Tuple{Float64, Float64}`: The cutoff radii for the first and second nearest neighbors.

# Returns
- `energies::Vector{Float64}`: A vector of the energies of each configuration.
- `configurations::Vector{Vector{Bool}}`: A vector of the configurations of the lattice system.

"""

function exact_enumeration(
    primitive_lattice_vectors::Matrix{Float64},
    basis::Vector{Tuple{Float64, Float64, Float64}},
    supercell_dimensions::Tuple{Int64, Int64, Int64},
    periodicity::Tuple{Bool, Bool, Bool},
    number_occupied_sites::Int64,
    adsorptions::Vector{Bool},
    adsorption_energy::Float64,
    nn_energy::Float64,
    nnn_energy::Float64,
    cutoff_radii::Tuple{Float64, Float64}
)
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

    # Generate LatticeSystem objects for each configuration
    lattices = [LatticeSystem(primitive_lattice_vectors, basis, supercell_dimensions, periodicity, occupations, adsorptions, cutoff_radii) for occupations in all_occupation_vectors]

    # Generate LatticeWalker objects for each lattice system
    walkers = [LatticeWalker(lattice) for lattice in lattices]

    # Compute energies for each walker
    for walker in walkers
        e_interaction = interaction_energy(walker.configuration, adsorption_energy, nn_energy, nnn_energy)
        walker.energy = e_interaction
    end

    # Extract energies and configurations
    energies = [walker.energy for walker in walkers]
    configurations = [walker.configuration for walker in walkers]

    return energies, configurations
end

"""
    nvt_monte_carlo(lattice::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, temperature::Float64, num_steps::Int64, random_seed::Int64)

Perform a Markov chain Monte Carlo simulation of a LatticeSystem in the NVT ensemble using the Metropolis-Hastings algorithm.

# Arguments
- `lattice::LatticeSystem`: The initial lattice configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.
- `temperature::Float64`: The temperature of the system.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `random_seed::Int64`: The seed for the random number generator.

# Returns
- `energies::Vector{Float64}`: A vector of the energies of the system at each step.
- `configurations::Vector{LatticeSystem}`: A vector of the configurations of the system at each step.
- `accepted_steps::Int64`: The number of accepted steps.
"""

function nvt_monte_carlo(
    lattice::LatticeSystem,
    adsorption_energy::Float64,
    nn_energy::Float64,
    nnn_energy::Float64,
    temperature::Float64,
    num_steps::Int64,
    random_seed::Int64
)
    # Set the random seed
    Random.seed!(random_seed)
    
    energies = Float64[]
    configurations = Vector{LatticeSystem}()
    accepted_steps = 0

    current_lattice = deepcopy(lattice)
    current_energy = interaction_energy(current_lattice, adsorption_energy, nn_energy, nnn_energy)
    
    push!(energies, current_energy)
    push!(configurations, deepcopy(current_lattice))

    for _ in 1:num_steps
        # Select a random site
        site_index = rand(1:length(current_lattice.occupations))
        
        # Propose a swap in occupation state (only if it maintains constant N)
        proposed_lattice = deepcopy(current_lattice)
        
        if proposed_lattice.occupations[site_index] == 1
            vacant_sites = findall(x -> x == 0, proposed_lattice.occupations)
            if length(vacant_sites) > 0
                swap_site = rand(vacant_sites)
                proposed_lattice.occupations[site_index] = 0
                proposed_lattice.occupations[swap_site] = 1
            else
                continue
            end
        else
            occupied_sites = findall(x -> x == 1, proposed_lattice.occupations)
            if length(occupied_sites) > 0
                swap_site = rand(occupied_sites)
                proposed_lattice.occupations[site_index] = 1
                proposed_lattice.occupations[swap_site] = 0
            else
                continue
            end
        end

        # Calculate the proposed energy
        proposed_energy = interaction_energy(proposed_lattice, adsorption_energy, nn_energy, nnn_energy)

        # Metropolis-Hastings acceptance criterion
        k_B = 8.617_333_262e-5  # eV K-1
        ΔE = proposed_energy - current_energy
        if ΔE < 0 || rand() < exp(-ΔE / (k_B * temperature))
            current_lattice = deepcopy(proposed_lattice)
            current_energy = proposed_energy
            accepted_steps += 1
        end

        push!(energies, current_energy)
        push!(configurations, deepcopy(current_lattice))
    end

    return energies, configurations, accepted_steps
end

"""
    wang_landau(lattice::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, 
                num_steps::Int64, flatness_criterion::Float64, f_initial::Float64, 
                f_min::Float64, energy_bins::Int64, random_seed::Int64)

Perform the Wang-Landau algorithm to compute the density of states for a LatticeSystem.

# Arguments
- `lattice::LatticeSystem`: The initial lattice configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `flatness_criterion::Float64`: The criterion for flatness of the histogram.
- `f_initial::Float64`: The initial modification factor.
- `f_min::Float64`: The minimum modification factor.
- `energy_bins::Int64`: The number of bins for the energy histogram.
- `random_seed::Int64`: The seed for the random number generator.

# Returns
- `density_of_states::Vector{Float64}`: The estimated density of states.
- `histogram::Vector{Int64}`: The histogram of visited energy states.
- `energies::Vector{Float64}`: The energies of the system at each step.
- `configurations::Vector{LatticeSystem}`: The configurations of the system at each step.
"""

function wang_landau(
    lattice::LatticeSystem,
    adsorption_energy::Float64,
    nn_energy::Float64,
    nnn_energy::Float64,
    num_steps::Int64,
    flatness_criterion::Float64,
    f_initial::Float64,
    f_min::Float64,
    energy_bins::Int64,
    random_seed::Int64
)
    # Set the random seed
    Random.seed!(random_seed)
    
    # Initialize density of states and histogram
    density_of_states = zeros(Float64, energy_bins)
    histogram = zeros(Int64, energy_bins)
    energies = Float64[]
    configurations = Vector{LatticeSystem}()
    
    f = f_initial
    current_lattice = deepcopy(lattice)
    current_energy = interaction_energy(current_lattice, adsorption_energy, nn_energy, nnn_energy)
    energy_min, energy_max = current_energy, current_energy

    push!(energies, current_energy)
    push!(configurations, deepcopy(current_lattice))

    # Define the binning function for energy
    bin_width = 1.0  # Initialize with a non-zero value to avoid division by zero
    get_bin = energy -> clamp(Int(floor((energy - energy_min) / bin_width)) + 1, 1, energy_bins)

    while f > f_min
        for i in 1:num_steps
            # Select a random site
            site_index = rand(1:length(current_lattice.occupations))
            
            # Propose a swap in occupation state (only if it maintains constant N)
            proposed_lattice = deepcopy(current_lattice)
            
            if proposed_lattice.occupations[site_index] == 1
                vacant_sites = findall(x -> x == 0, proposed_lattice.occupations)
                if length(vacant_sites) > 0
                    swap_site = rand(vacant_sites)
                    proposed_lattice.occupations[site_index] = 0
                    proposed_lattice.occupations[swap_site] = 1
                else
                    continue
                end
            else
                occupied_sites = findall(x -> x == 1, proposed_lattice.occupations)
                if length(occupied_sites) > 0
                    swap_site = rand(occupied_sites)
                    proposed_lattice.occupations[site_index] = 1
                    proposed_lattice.occupations[swap_site] = 0
                else
                    continue
                end
            end

            # Calculate the proposed energy
            proposed_energy = interaction_energy(proposed_lattice, adsorption_energy, nn_energy, nnn_energy)
            
            # Ensure the energy is within bounds
            energy_min = min(energy_min, proposed_energy)
            energy_max = max(energy_max, proposed_energy)
            bin_width = (energy_max - energy_min) / energy_bins

            current_bin = get_bin(current_energy)
            proposed_bin = get_bin(proposed_energy)

            # Wang-Landau acceptance criterion
            if density_of_states[proposed_bin] == 0 || rand() < exp(density_of_states[current_bin] - density_of_states[proposed_bin])
                current_lattice = deepcopy(proposed_lattice)
                current_energy = proposed_energy
            end

            current_bin = get_bin(current_energy)
            density_of_states[current_bin] += log(f)
            histogram[current_bin] += 1

            push!(energies, current_energy)
            push!(configurations, deepcopy(current_lattice))
        end

        # Check for flatness of the histogram
        if minimum(histogram) > flatness_criterion * mean(histogram)
            f = sqrt(f)
            histogram .= 0
        end
    end

    return density_of_states, histogram, energies, configurations
end