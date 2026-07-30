"""
    _check_shell_counts(lattice_neighbors, n_coupled::Int)

Throw an `ArgumentError` when a Hamiltonian couples more neighbor shells
than the lattice's neighbor list provides. Called from the energy kernels,
so it is a single integer comparison in the common case.
"""
@inline function _check_shell_counts(lattice_neighbors::Vector{Vector{Vector{Int64}}}, n_coupled::Int)
    n_shells = isempty(lattice_neighbors) ? 0 : length(lattice_neighbors[1])
    if n_coupled > n_shells
        throw(ArgumentError(
            "the Hamiltonian couples $n_coupled neighbor shells but the " *
            "lattice provides only $n_shells (= length(cutoff_radii)); " *
            "extend cutoff_radii so every coupled shell exists"))
    end
    return nothing
end

"""
    lattice_interaction_energy(lattice_occupations::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U})

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice_occupations::Vector{Bool}`: The lattice occupation configuration.
- `lattice_neighbors::Vector{Vector{Vector{Int64}}}`: The lattice neighbor list.
- `h::GenericLatticeHamiltonian{N,U}`: The generic lattice Hamiltonian parameters.

# Returns
- `e_interaction::U`: The interaction energy of the lattice configuration.

Throws an `ArgumentError` when the Hamiltonian couples more neighbor shells
than the lattice's neighbor list provides (`N > length(cutoff_radii)`),
which would otherwise surface as a raw `BoundsError` from the innermost
loop. The converse mismatch (fewer coupled shells than the lattice carries)
is legal — the outer shells are simply not coupled — and is flagged once,
with a warning, at `LatticeGasWalkers` construction.
"""
function lattice_interaction_energy(lattice_occupations::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U}) where {N,U}
    _check_shell_counts(lattice_neighbors, N)
    e_interaction::U = 0.0*unit(h.on_site_interaction)
    for index in eachindex(lattice_occupations)
        if lattice_occupations[index]
            for n in 1:N
                for neighbor in lattice_neighbors[index][n]
                    if lattice_occupations[neighbor]
                        e_interaction += h.nth_neighbor_interactions[n] / 2
                    end
                end
            end
        end
    end
    return e_interaction
end

"""
    inter_component_energy(lattice1::Vector{Bool}, lattice2::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U})

Compute the interaction energy between two lattice configurations using the Hamiltonian parameters.

# Arguments
- `lattice1::Vector{Bool}`: The first lattice configuration.
- `lattice2::Vector{Bool}`: The second lattice configuration.
- `lattice_neighbors::Vector{Vector{Vector{Int64}}}`: The lattice neighbor list.
- `h::GenericLatticeHamiltonian{N,U}`: The generic lattice Hamiltonian parameters.

# Returns
- `e_interaction::U`: The interaction energy between the two lattice configurations.

Throws an `ArgumentError` on the same shell-count mismatch as
[`lattice_interaction_energy`](@ref).
"""
function inter_component_energy(lattice1::Vector{Bool}, lattice2::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U}) where {N,U}
    _check_shell_counts(lattice_neighbors, N)
    e_interaction::U = 0.0*unit(h.on_site_interaction)
    for index in eachindex(lattice1)
        if lattice1[index]
            for n in 1:N
                for neighbor in lattice_neighbors[index][n]
                    if lattice2[neighbor]
                        e_interaction += h.nth_neighbor_interactions[n]
                    end
                end
            end
        end
    end
    return e_interaction
end

"""
    interacting_energy(lattice::SLattice, h::GenericLatticeHamiltonian{N})

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice::SLattice`: The lattice configuration.
- `h::GenericLatticeHamiltonian{N,U}`: The generic lattice Hamiltonian parameters.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""
function interacting_energy(lattice::SLattice, h::GenericLatticeHamiltonian{N,U}) where {N,U}
    # for SLattice with a single-component Hamiltonian
    e_interaction::U = lattice_interaction_energy(lattice.components[1], lattice.neighbors, h)
    e_adsorption::U = sum(lattice.components[1] .& lattice.adsorptions) * h.on_site_interaction
    return e_interaction + e_adsorption
end

function interacting_energy(lattice::SLattice, h::MLatticeHamiltonian{C,N,U}) where {C,N,U}
    # for SLattice with a multi-component Hamiltonian, taking the first element of the Hamiltonian matrix
    ham = h.Hamiltonians[1,1]
    e_interaction::U = lattice_interaction_energy(lattice.components[1], lattice.neighbors, ham)
    e_adsorption::U = sum(lattice.components[1] .& lattice.adsorptions) * ham.on_site_interaction
    return e_interaction + e_adsorption
end

"""
    interacting_energy(lattice::MLattice{C,G}, h::MLatticeHamiltonian{C,N,U})

Compute the interaction energy of a multi-component lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice::MLattice{C,G}`: The multi-component lattice configuration.
- `h::MLatticeHamiltonian{C,N,U}`: The multi-component lattice Hamiltonian parameters.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""
function interacting_energy(lattice::MLattice{C,G}, h::MLatticeHamiltonian{C,N,U}) where {C,G,N,U}
    # for MLattice with a multi-component Hamiltonian
    adsorption_energy = 0.0*unit(h.Hamiltonians[1,1].on_site_interaction)
    interaction_energy = 0.0*unit(h.Hamiltonians[1,1].on_site_interaction)
    for i in 1:C
        interaction_energy += lattice_interaction_energy(lattice.components[i], lattice.neighbors, h.Hamiltonians[i,i])
        adsorption_energy += sum(lattice.components[i] .& lattice.adsorptions) * h.Hamiltonians[i,i].on_site_interaction
        for j in (i+1):C
            interaction_energy += inter_component_energy(lattice.components[i], lattice.components[j], lattice.neighbors, h.Hamiltonians[i,j])
        end
    end
    @debug "interaction energy: $interaction_energy, adsorption energy: $adsorption_energy"
    return interaction_energy + adsorption_energy
end