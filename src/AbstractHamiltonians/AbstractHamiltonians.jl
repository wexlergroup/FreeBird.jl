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
    println(io, "    nth_neighbor_interactions: ", [hamiltonian.nth_neighbor_interactions[i].val for i in 1:N], " ", unit(U))
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

end # module Hamiltonians