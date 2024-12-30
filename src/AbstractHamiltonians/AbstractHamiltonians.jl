"""
Method for defining and implementing Hamiltonians.
"""
module AbstractHamiltonians

using Unitful
using StaticArrays
using StaticArrays

export ClassicalHamiltonian
export LatticeGasHamiltonian, GenericLatticeHamiltonian, MLatticeHamiltonian

abstract type ClassicalHamiltonian end

"""
    struct LatticeGasHamiltonian <: ClassicalHamiltonian

The `LatticeGasHamiltonian` struct represents the parameters for a lattice-gas Hamiltonian.

# Fields
- `adsorption_energy::typeof(1.0u"eV")`: The energy of adsorption on the lattice.
- `nn_interaction_energy::typeof(1.0u"eV")`: The energy of nearest-neighbor interactions.
- `nnn_interaction_energy::typeof(1.0u"eV")`: The energy of next-nearest-neighbor interactions.

"""

struct LatticeGasHamiltonian <: ClassicalHamiltonian
    adsorption_energy::typeof(1.0u"eV")
    nn_interaction_energy::typeof(1.0u"eV")
    nnn_interaction_energy::typeof(1.0u"eV")
end

function Base.show(io::IO, hamiltonian::LatticeGasHamiltonian)
    println(io, "LatticeGasHamiltonian:")
    println(io, "    adsorption_energy:      ", hamiltonian.adsorption_energy)
    println(io, "    nn_interaction_energy:  ", hamiltonian.nn_interaction_energy)
    println(io, "    nnn_interaction_energy: ", hamiltonian.nnn_interaction_energy)
end


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
julia> hams = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4]
julia> mlham = MLatticeHamiltonian(hams)
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
        return MLatticeHamiltonian{c,N,U}(reshape(hams, c, c))
    elseif length(hams) == c*(c+1)/2
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