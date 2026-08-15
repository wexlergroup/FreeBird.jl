"""
    AbstractHamiltonians

Module for defining and implementing Hamiltonians.
"""
module AbstractHamiltonians

using Unitful
using StaticArrays
using StaticArrays

export AbstractHamiltonian
export ClassicalHamiltonian
export GenericLatticeHamiltonian, MLatticeHamiltonian
export ClusterInteraction, ClusterLatticeHamiltonian
export SiteFieldLatticeHamiltonian
export supports_site_deltas

abstract type AbstractHamiltonian end

abstract type ClassicalHamiltonian <: AbstractHamiltonian end


"""
    struct GenericLatticeHamiltonian{N,U} <: ClassicalHamiltonian

The `GenericLatticeHamiltonian` struct represents a generic lattice Hamiltonian. 
It has an on-site interaction energy and a `N`-elements vector of nth-neighbor interaction energies.
Units of energy `U` is also specified.

# Fields
- `on_site_interaction::U`: The energy of on-site interactions.
- `nth_neighbor_interactions::SVector{N, U}`: The energy of nth-neighbor interactions.

# Constructors
```julia
GenericLatticeHamiltonian(on_site_interaction::Float64, nth_neighbor_interactions::Vector{Float64}, energy_units::Unitful.Units)
GenericLatticeHamiltonian(on_site_interaction::U, nth_neighbor_interactions::Vector{U}) where U
```

All couplings must be finite; `Inf` is rejected by the constructor because it
stalls nested sampling silently (every ceiling comparison becomes `Inf >= Inf`).

## Hard-core (athermal) lattice models

Nearest-neighbor exclusion models — hard squares on the square lattice, hard
hexagons on the triangular lattice — are expressed with a *finite* repulsive
coupling on a lattice built with a single cutoff shell, e.g.
`GenericLatticeHamiltonian(0.0, [1.0], u"eV")` with `cutoff_radii=[1.1]`.
The energy is then `J × (number of excluded-neighbor pairs)`, an
integer-leveled ladder whose `E = 0` manifold is exactly the hard-core
configuration space, and the nested-sampling descent through the violating
levels measures the hard-core partition function against the full,
unrestricted prior (the `(1+z0)^M` normalization stays valid). Two numerical
windows apply:

- **Sampling**: choose the NS `energy_perturbation` δ with
  `eps(E_max) ≪ δ ≪ J`, where `E_max = J·M·c/2` is the maximal violation
  energy on `M` sites with excluded-shell coordination `c`. Concretely, keep
  `δ / eps(E_max)` above ~`K²` so the `K` walkers draw distinct tie-breaking
  values on a plateau: δ = `1e-9` at `J = 1 eV` is safe through
  `E_max ≈ 10³ eV` (hard squares near `M ≈ 500`, hard hexagons near
  `M ≈ 330`); scale δ proportionally for larger lattices.
- **Post-processing**: evaluate at a temperature low enough that
  `β·J ≥ (ladder depth in nats) + ~40`, with the depth
  `n_iters · ln((K + n_cull)/K)`, while keeping `β·δ ≪ 1`; every violating
  level's Boltzmann factor is then an exact zero to double precision, the
  allowed levels' deviate from 1 only by `O(β·δ)`, and all observables are
  athermal (functions of the activity `z = exp(βμ)` alone).

## Examples
```jldoctest
julia> ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
GenericLatticeHamiltonian{2,Quantity{Float64, 𝐋² 𝐌 𝐓⁻², Unitful.FreeUnits{(eV,), 𝐋² 𝐌 𝐓⁻², nothing}}}:
    on_site_interaction:      -0.04 eV
    nth_neighbor_interactions: [-0.01, -0.0025] eV


julia> ham = GenericLatticeHamiltonian(-0.04u"eV", [-0.01, -0.0025]*u"eV")
GenericLatticeHamiltonian{2,Quantity{Float64, 𝐋² 𝐌 𝐓⁻², Unitful.FreeUnits{(eV,), 𝐋² 𝐌 𝐓⁻², nothing}}}:
    on_site_interaction:      -0.04 eV
    nth_neighbor_interactions: [-0.01, -0.0025] eV
```
"""
struct GenericLatticeHamiltonian{N,U} <: ClassicalHamiltonian
    on_site_interaction::U
    nth_neighbor_interactions::SVector{N, U}
    function GenericLatticeHamiltonian{N,U}(on_site_interaction::U, nth_neighbor_interactions::Vector{U}) where {N,U}
        if length(nth_neighbor_interactions) != N
            throw(ArgumentError("Length of nth_neighbor_interactions must match the number of neighbors"))
        end
        if !isfinite(on_site_interaction) || any(!isfinite, nth_neighbor_interactions)
            throw(ArgumentError("non-finite couplings are not supported: an Inf coupling " *
                "makes every nested-sampling ceiling comparison degenerate (Inf >= Inf) and " *
                "the sampler stalls silently. Model hard-core exclusion with a finite " *
                "repulsive coupling instead — see the hard-core recipe in the " *
                "GenericLatticeHamiltonian docstring."))
        end
        return new(on_site_interaction, SVector{N, U}(nth_neighbor_interactions))
    end
end

function GenericLatticeHamiltonian(on_site_interaction::U, nth_neighbor_interactions::Vector{U}) where U
    return GenericLatticeHamiltonian{length(nth_neighbor_interactions),U}(on_site_interaction, nth_neighbor_interactions)
end

function GenericLatticeHamiltonian(on_site_interaction::Float64, nth_neighbor_interactions::Vector{Float64}, energy_units::Unitful.Units)
    return GenericLatticeHamiltonian{length(nth_neighbor_interactions),typeof(0.0*energy_units)}(on_site_interaction * energy_units, nth_neighbor_interactions * energy_units)
end

function Base.show(io::IO, hamiltonian::GenericLatticeHamiltonian{N,U}) where {N,U}
    println(io, "GenericLatticeHamiltonian{$N,$U}:")
    println(io, "    on_site_interaction:      ", hamiltonian.on_site_interaction)
    println(io, "    nth_neighbor_interactions: ", [ustrip(hamiltonian.nth_neighbor_interactions[i]) for i in 1:N], " ", unit(U))
end
    
    """
        struct MLatticeHamiltonian{C,N,U} <: ClassicalHamiltonian

The `MLatticeHamiltonian` struct represents a multi-component lattice Hamiltonian.
It has a matrix of `GenericLatticeHamiltonian{N,U}`.

# Fields
- `Hamiltonians::Matrix{GenericLatticeHamiltonian{N,U}}`: The matrix of `GenericLatticeHamiltonian{N,U}`.

# Constructors
```julia
MLatticeHamiltonian(Hamiltonians::Vector{GenericLatticeHamiltonian{N,U}})
```
## Examples
```jldoctest
julia> hams = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4] # full flattened matrix
julia> mlham = MLatticeHamiltonian(hams)
julia> hams = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:3] # upper triangular matrix
julia> mlham = MLatticeHamiltonian(2, hams)
```    
"""
struct MLatticeHamiltonian{C,N,U} <: ClassicalHamiltonian
    Hamiltonians::Matrix{GenericLatticeHamiltonian{N,U}}

    function MLatticeHamiltonian{C,N,U}(Hamiltonians::Matrix{GenericLatticeHamiltonian{N,U}}) where {C,N,U}
        if size(Hamiltonians) != (C, C)
            throw(ArgumentError("the size of the matrix is not compatible with the number of components."))
        end
        return new(Hamiltonians)
    end
end

function MLatticeHamiltonian(c::Int, hams::Vector{GenericLatticeHamiltonian{N,U}}) where {N,U}
    if length(hams) == c^2
        @info "Creating MLatticeHamiltonian from a flattened matrix. By specifying
        $(length(hams)) sets of GenericLatticeHamiltonian, a $c x $c matrix is constructed. If this 
        was not your intention, please check the documentation or raise an issue."
        return MLatticeHamiltonian{c,N,U}(reshape(hams, c, c))
    elseif length(hams) == c*(c+1)/2
        @info "Creating MLatticeHamiltonian from a vector of GenericLatticeHamiltonian. By specifying
        $(length(hams)) sets of GenericLatticeHamiltonian, a $c x $c matrix is constructed. If this
        was not your intention, please check the documentation or raise an issue."
        ham_matrix = Matrix{GenericLatticeHamiltonian{N,U}}(undef, c, c)
        k = 1
        for i in 1:c
            for j in i:c
                ham_matrix[i,j] = hams[k]
                ham_matrix[j,i] = hams[k]
                k += 1
            end
        end
        return MLatticeHamiltonian{c,N,U}(ham_matrix)
    else
        throw(ArgumentError("the number of elements in the vector must be equal to c^2 or c*(c+1)/2."))
    end
    
end

function Base.show(io::IO, hamiltonian::MLatticeHamiltonian{C,N,U}) where {C,N,U}
    println(io, "MLatticeHamiltonian{$C,$N,$U}:")
    for i in 1:C
       for j in 1:C
           println(io, "    Hamiltonians[$i, $j]: ", hamiltonian.Hamiltonians[i,j])
       end
    end
end

"""
    struct ClusterInteraction{K,U}

A single multi-body (cluster) interaction figure on a lattice: one coupling
plus the explicit list of the figure's embeddings, each an ordered `K`-tuple
of site indices. The coupling is the energy contributed per fully occupied
embedding. Each unordered site set appears exactly once, in canonical
(strictly increasing) form, so double counting of symmetric motifs is
structurally impossible. Embedding lists are typically produced by
`enumerate_motif_embeddings` under the torus (minimum-image) convention,
matching how the pair shells count neighbors.

Only orders `K ≥ 3` are represented: on-site and pair terms belong in the
wrapped `GenericLatticeHamiltonian` of a [`ClusterLatticeHamiltonian`](@ref).

Note: "cluster" here names an interaction figure of a cluster expansion,
not the geometric cluster Monte Carlo moves elsewhere in FreeBird.

# Constructor
```julia
ClusterInteraction(coupling, embeddings::Vector{NTuple{K,Int}})
```
Throws `ArgumentError` on `K < 3`, a non-finite coupling, an embedding that
is not strictly increasing or contains a non-positive index, or a duplicate
embedding; warns on an empty embedding list (which contributes exactly
zero energy).
"""
struct ClusterInteraction{K,U}
    coupling::U
    embeddings::Vector{NTuple{K,Int}}

    function ClusterInteraction{K,U}(coupling::U, embeddings::Vector{NTuple{K,Int}}) where {K,U}
        if K < 3
            throw(ArgumentError(
                "ClusterInteraction is for 3-body and higher motifs (got K = $K); " *
                "on-site and pair terms belong in the wrapped GenericLatticeHamiltonian"))
        end
        if !isfinite(coupling)
            throw(ArgumentError(
                "non-finite cluster couplings are not supported; model hard " *
                "constraints with finite couplings (see the hard-core recipe " *
                "in the GenericLatticeHamiltonian docstring)"))
        end
        isempty(embeddings) &&
            @warn "ClusterInteraction with an empty embedding list contributes exactly zero energy"
        seen = Set{NTuple{K,Int}}()
        for e in embeddings
            if e[1] < 1
                throw(ArgumentError("embedding $e contains a non-positive site index"))
            end
            for k in 1:K-1
                if e[k] >= e[k+1]
                    throw(ArgumentError(
                        "embedding $e is not in canonical (strictly increasing) " *
                        "order; each unordered site set must appear exactly once"))
                end
            end
            if e in seen
                throw(ArgumentError("duplicate embedding $e would double-count its energy"))
            end
            push!(seen, e)
        end
        return new(coupling, embeddings)
    end
end

function ClusterInteraction(coupling::U, embeddings::Vector{NTuple{K,Int}}) where {K,U}
    return ClusterInteraction{K,U}(coupling, embeddings)
end

function Base.show(io::IO, c::ClusterInteraction{K,U}) where {K,U}
    print(io, "ClusterInteraction{$K}: coupling = ", c.coupling,
          ", ", length(c.embeddings), " embeddings")
end

"""
    struct ClusterLatticeHamiltonian{N,U} <: ClassicalHamiltonian

A lattice Hamiltonian with on-site, pair-shell, and multi-body (cluster)
terms: a `GenericLatticeHamiltonian{N,U}` carrying the on-site energy and
the `N` pair-shell couplings, plus a vector of [`ClusterInteraction`](@ref)
figures of arbitrary orders `K ≥ 3`. The total energy is the pair part
plus, for every figure, the coupling times the number of its fully
occupied embeddings.

Sign and unit conventions follow `GenericLatticeHamiltonian`: positive
couplings are repulsive. Published binding-energy parameter sets (e.g. the
O-Pd(100) lattice-gas Hamiltonian of Zhang, Blum & Reuter, PRB 75, 235406
(2007), whose optimum expansion carries four pair shells, four trios, and
one quattro) map onto FreeBird by negating every parameter.

Evaluated by `interacting_energy(lattice, h)` for single-component
(`SLattice`) configurations, and therefore usable in every sampler that
calls it generically.

# Constructor
```julia
ClusterLatticeHamiltonian(pair_ham::GenericLatticeHamiltonian{N,U}, clusters)
```
`clusters` is a vector of `ClusterInteraction`s whose coupling type matches
`U`; heterogeneous orders `K` are allowed.
"""
struct ClusterLatticeHamiltonian{N,U} <: ClassicalHamiltonian
    pair_ham::GenericLatticeHamiltonian{N,U}
    clusters::Vector{ClusterInteraction{K,U} where K}

    function ClusterLatticeHamiltonian(pair_ham::GenericLatticeHamiltonian{N,U},
                                       clusters::AbstractVector) where {N,U}
        converted = Vector{ClusterInteraction{K,U} where K}()
        for (i, c) in enumerate(clusters)
            if !(c isa ClusterInteraction)
                throw(ArgumentError(
                    "clusters[$i] is a $(typeof(c)); expected a ClusterInteraction"))
            end
            if !(typeof(c.coupling) == U)
                throw(ArgumentError(
                    "clusters[$i] has coupling type $(typeof(c.coupling)), which " *
                    "does not match the pair Hamiltonian's $U; construct every " *
                    "coupling with the same units and value type"))
            end
            push!(converted, c)
        end
        return new{N,U}(pair_ham, converted)
    end
end

function Base.show(io::IO, h::ClusterLatticeHamiltonian{N,U}) where {N,U}
    println(io, "ClusterLatticeHamiltonian{$N,$U}:")
    println(io, "    pair part: ", h.pair_ham)
    println(io, "    clusters:")
    for c in h.clusters
        println(io, "        ", c)
    end
end

"""
    _coupling_type(hamiltonian)

Internal trait returning the concrete coupling type `U` of a lattice
Hamiltonian (the exact `typeof` shared by its on-site, pair, and cluster
couplings), or `nothing` for Hamiltonian types that do not declare one. The
[`SiteFieldLatticeHamiltonian`](@ref) constructor uses it to require that
every field entry's type equals the base's coupling type (same units and
value type; no promotion), mirroring the per-cluster coupling check of
[`ClusterLatticeHamiltonian`](@ref). New `ClassicalHamiltonian` types opt
in by adding a method.
"""
_coupling_type(hamiltonian) = nothing
_coupling_type(::GenericLatticeHamiltonian{N,U}) where {N,U} = U
_coupling_type(::MLatticeHamiltonian{C,N,U}) where {C,N,U} = U
_coupling_type(::ClusterLatticeHamiltonian{N,U}) where {N,U} = U

"""
    struct SiteFieldLatticeHamiltonian{H,U} <: ClassicalHamiltonian

A wrapper Hamiltonian adding a site-resolved on-site energy (a "field") to
a lattice Hamiltonian: the total energy is the wrapped base Hamiltonian's
energy plus `Σ field[i]` over every occupied site `i`. One number per site
breaks the site equivalence the base Hamiltonians assume, which is what
layered adsorption physics needs: a substrate potential decaying with
height, canonically a contact term plus an inverse-cube tail in the
multilayer lattice-gas models of de Oliveira and Griffiths
[Surf. Sci. 71, 687 (1978)] and Pandit, Schick and Wortis
[Phys. Rev. B 26, 5112 (1982)], is expressed as a per-layer profile
broadcast over sites with `layer_field`.

# Fields
- `base::H`: the wrapped Hamiltonian (`GenericLatticeHamiltonian`,
  `MLatticeHamiltonian`, or `ClusterLatticeHamiltonian`).
- `field::Vector{U}`: one on-site energy per site, added once per occupied
  site; `U` equals the base's coupling type exactly.

# Additive contract

The field adds on top of the base's own on-site term: a base with
`on_site_interaction = ε` still contributes `ε × (number of occupied
adsorption sites)` through `lattice.adsorptions`. The two on-site channels
compose additively, so supplying a full per-site field while the base
carries a nonzero on-site energy double-counts every occupied adsorption
site. When a full field is supplied, zero the base's on-site energy and
put the entire on-site physics in the field.

The wrapper subsumes the adsorption mechanism exactly: for a lattice with
adsorption mask `A = lattice.adsorptions`,

    GenericLatticeHamiltonian(ε, J)   # ε once per occupied adsorption site
    ≡ SiteFieldLatticeHamiltonian(GenericLatticeHamiltonian(zero(ε), J), ε .* A)

since the masked field `ε .* A` carries `ε` on each adsorption site and an
exact zero elsewhere.

## Layer-profile recipe

    lat = SLattice{SquareLattice}(supercell_dimensions=(4, 4, 3),
                                  periodicity=(true, true, false),
                                  cutoff_radii=[1.1])
    field = layer_field(lat, [-0.27, -0.03375, -0.01] .* u"eV")  # inverse-cube tail
    h = SiteFieldLatticeHamiltonian(GenericLatticeHamiltonian(0.0, [-0.01], u"eV"), field)

# Constructors
```julia
SiteFieldLatticeHamiltonian(base::ClassicalHamiltonian, field::Vector)
SiteFieldLatticeHamiltonian(base::ClassicalHamiltonian, field::Vector{Float64}, energy_units::Unitful.Units)
```

The three-argument form attaches `energy_units` to a plain `Float64`
vector, like the units convenience constructor of
[`GenericLatticeHamiltonian`](@ref). The field vector is copied.

Throws an `ArgumentError` when the base is itself a
`SiteFieldLatticeHamiltonian` (compose site fields by adding the vectors,
not by nesting wrappers), when the field is empty, when the base's
coupling type cannot be determined, when any entry's type differs from the
base's coupling type (exact `typeof` equality: same units and value type,
no promotion), or when any entry is non-finite (an `Inf` makes every
nested-sampling ceiling comparison degenerate and the sampler stalls
silently; see the hard-core recipe in the
[`GenericLatticeHamiltonian`](@ref) docstring).

The field length is *not* checked here (no lattice is in scope at
construction) but once per run at setup, in the `LatticeGasWalkers`
constructor and the raw-lattice `wang_landau`/`nvt_monte_carlo` entry
points, where a length differing from the lattice's site count throws.
"""
struct SiteFieldLatticeHamiltonian{H<:ClassicalHamiltonian,U} <: ClassicalHamiltonian
    base::H
    field::Vector{U}

    function SiteFieldLatticeHamiltonian(base::ClassicalHamiltonian, field::AbstractVector)
        if base isa SiteFieldLatticeHamiltonian
            throw(ArgumentError(
                "nesting SiteFieldLatticeHamiltonian wrappers is not supported; " *
                "compose site fields by adding the vectors: " *
                "SiteFieldLatticeHamiltonian(base.base, base.field .+ field)"))
        end
        if isempty(field)
            throw(ArgumentError(
                "field is empty; supply one on-site energy per lattice site " *
                "(build a layer profile with layer_field)"))
        end
        U = _coupling_type(base)
        if U === nothing
            throw(ArgumentError(
                "cannot determine the coupling type of a $(typeof(base)); " *
                "SiteFieldLatticeHamiltonian supports GenericLatticeHamiltonian, " *
                "MLatticeHamiltonian and ClusterLatticeHamiltonian bases, or any " *
                "Hamiltonian type that defines AbstractHamiltonians._coupling_type"))
        end
        for (i, v) in enumerate(field)
            if !(typeof(v) == U)
                throw(ArgumentError(
                    "field[$i] has type $(typeof(v)), which does not match the " *
                    "base Hamiltonian's coupling type $U; construct every field " *
                    "entry with the same units and value type as the base couplings"))
            end
            if !isfinite(v)
                throw(ArgumentError(
                    "field[$i] is not finite; non-finite on-site energies are not " *
                    "supported: an Inf makes every nested-sampling ceiling " *
                    "comparison degenerate (Inf >= Inf) and the sampler stalls " *
                    "silently. Model site exclusion with a finite repulsive value " *
                    "instead; see the hard-core recipe in the " *
                    "GenericLatticeHamiltonian docstring."))
            end
        end
        return new{typeof(base),U}(base, collect(U, field))
    end
end

function SiteFieldLatticeHamiltonian(base::ClassicalHamiltonian, field::Vector{Float64}, energy_units::Unitful.Units)
    return SiteFieldLatticeHamiltonian(base, field .* energy_units)
end

_coupling_type(::SiteFieldLatticeHamiltonian{H,U}) where {H,U} = U

function Base.show(io::IO, h::SiteFieldLatticeHamiltonian{H,U}) where {H,U}
    println(io, "SiteFieldLatticeHamiltonian{$H,$U}:")
    lo, hi = extrema(h.field)
    println(io, "    field: ", length(h.field), " sites, extrema [",
            ustrip(lo), ", ", ustrip(hi), "] ", unit(U))
    println(io, "    base: ", h.base)
end

"""
    supports_site_deltas(hamiltonian::ClassicalHamiltonian) -> Bool

Opt-in trait: whether an exact O(z) single-site occupancy-flip energy delta
(`site_flip_delta` in `EnergyEval`) is available for the Hamiltonian type.

The `false` fallback is a contract, not a failure: consumers fall back to
full recomputation instead of silently mis-evaluating. `ClusterLatticeHamiltonian`
deliberately stays `false` (an exact cluster delta needs precomputed
site-to-embedding incidence lists that do not exist yet), and the site-field
wrapper delegates to its base. New `ClassicalHamiltonian` types opt in by
adding a method alongside their `site_flip_delta` method.
"""
supports_site_deltas(::ClassicalHamiltonian) = false
supports_site_deltas(::GenericLatticeHamiltonian) = true
supports_site_deltas(::MLatticeHamiltonian) = true
supports_site_deltas(h::SiteFieldLatticeHamiltonian) = supports_site_deltas(h.base)

end # module Hamiltonians