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
    cluster_energy(occupations::Vector{Bool}, c::ClusterInteraction{K,U})

Energy contribution of one cluster figure: the coupling times the number of
its embeddings whose `K` sites are all occupied. The per-figure function
barrier keeps the inner loop type-stable across the heterogeneous cluster
orders of a `ClusterLatticeHamiltonian`. Bounds checks stay on: an
embedding referencing a site beyond the occupation vector must raise a
`BoundsError` rather than read memory (the liveset constructor validates
this once up front).
"""
function cluster_energy(occupations::Vector{Bool}, c::ClusterInteraction{K,U}) where {K,U}
    n = 0
    for emb in c.embeddings
        occupied = true
        for s in emb
            if !occupations[s]
                occupied = false
                break
            end
        end
        n += occupied ? 1 : 0
    end
    return n * c.coupling
end

"""
    interacting_energy(lattice::SLattice, h::ClusterLatticeHamiltonian{N,U})

Total energy under a multi-body lattice Hamiltonian: the wrapped pair
part (on-site + `N` pair shells, evaluated exactly as for a bare
`GenericLatticeHamiltonian`) plus every cluster figure's contribution.
Single-component (`SLattice`) configurations only.
"""
function interacting_energy(lattice::SLattice, h::ClusterLatticeHamiltonian{N,U}) where {N,U}
    e = interacting_energy(lattice, h.pair_ham)
    occ = lattice.components[1]
    for c in h.clusters
        e += cluster_energy(occ, c)
    end
    return e
end

"""
    site_field_energy(occupations::Vector{Bool}, field::Vector{U})

Energy contribution of a per-site field: the sum of `field[i]` over every
occupied site `i`. Iteration runs over `eachindex(occupations, field)`, so
a field whose length differs from the occupation vector (a Hamiltonian
built for a different lattice) raises a `DimensionMismatch` rather than
silently summing a truncated or padded range; the liveset constructor and
the raw-lattice sampler entry points validate the length once up front,
with a descriptive error.
"""
function site_field_energy(occupations::Vector{Bool}, field::Vector{U}) where U
    e::U = zero(U)
    for i in eachindex(occupations, field)
        if occupations[i]
            e += field[i]
        end
    end
    return e
end

"""
    interacting_energy(lattice::SLattice, h::SiteFieldLatticeHamiltonian{H,U})

Total energy under a site-field wrapper: the wrapped base Hamiltonian's
energy, evaluated exactly as if the base were passed directly (including
its `on_site_interaction × occupied-adsorption-sites` term), plus the
occupation-masked field sum of [`site_field_energy`](@ref).
Single-component (`SLattice`) configurations only.
"""
function interacting_energy(lattice::SLattice, h::SiteFieldLatticeHamiltonian{H,U}) where {H,U}
    e_base::U = interacting_energy(lattice, h.base)
    e_field::U = site_field_energy(lattice.components[1], h.field)
    return e_base + e_field
end

"""
    site_flip_delta(lattice::SLattice, h::GenericLatticeHamiltonian{N,U}, site::Int) where {N,U}

Exact energy change from flipping the occupancy of `site`:
`E(after) − E(before)`, computed as a local O(z) sum over the flipped site's
neighbor entries (z the per-site neighbor entry count), the lattice
counterpart of the audited `single_site_energy` path on the continuous side.

The formula respects both conventions of the full sweep
([`lattice_interaction_energy`](@ref)): each ordered neighbor entry carries
half a coupling, so an ordinary neighbor `j != site` (appearing once in the
flipped site's list and once in `j`'s) contributes `occ[j]` times the full
coupling to the delta, while a self-image entry (`j == site`, present under
`image_multiplicity = true`) appears only in the flipped site's own row and
contributes exactly half the coupling per image. The adsorption term mirrors
the on-site channel of the full evaluation.

Availability for a Hamiltonian type is declared by the
`supports_site_deltas` trait; consumers fall back to full recomputation when
it is `false`. Throws the same shell-count `ArgumentError` as the full sweep
when the Hamiltonian couples more shells than the lattice provides.
"""
function site_flip_delta(lattice::SLattice, h::GenericLatticeHamiltonian{N,U},
                         site::Int) where {N,U}
    _check_shell_counts(lattice.neighbors, N)
    occ = lattice.components[1]
    sgn = occ[site] ? -1 : 1
    acc::U = lattice.adsorptions[site] ? h.on_site_interaction :
             zero(h.on_site_interaction)
    nbrs = lattice.neighbors[site]
    for n in 1:N
        Jn = h.nth_neighbor_interactions[n]
        for j in nbrs[n]
            if j == site
                acc += Jn / 2
            elseif occ[j]
                acc += Jn
            end
        end
    end
    return sgn * acc
end

"""
    site_flip_delta(lattice::SLattice, h::MLatticeHamiltonian{C,N,U}, site::Int) where {C,N,U}

Single-component delta under a multi-component Hamiltonian: delegates to
`h.Hamiltonians[1, 1]`, matching the corresponding `interacting_energy`
method for `SLattice`.
"""
function site_flip_delta(lattice::SLattice, h::MLatticeHamiltonian{C,N,U},
                         site::Int) where {C,N,U}
    return site_flip_delta(lattice, h.Hamiltonians[1, 1], site)
end

"""
    site_flip_delta(lattice::SLattice, h::SiteFieldLatticeHamiltonian{H,U}, site::Int) where {H,U}

Delta under the site-field wrapper: the base Hamiltonian's delta plus the
signed field entry of the flipped site, mirroring the wrapper's
`interacting_energy` method.
"""
function site_flip_delta(lattice::SLattice, h::SiteFieldLatticeHamiltonian{H,U},
                         site::Int) where {H,U}
    base_delta = site_flip_delta(lattice, h.base, site)
    sgn = lattice.components[1][site] ? -1 : 1
    return base_delta + sgn * h.field[site]
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