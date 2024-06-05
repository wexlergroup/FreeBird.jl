"""
Module for defining and implementing potentials.
"""
module Potentials

using Unitful

export LJParameters, lj_energy
export CompositeLJParameters

abstract type LennardJonesParametersSets end


"""
    struct LJParameters

The `LJParameters` struct represents the parameters for the Lennard-Jones potential.

# Fields
- `epsilon::typeof(1.0u"eV")`: The energy scale of the potential.
- `sigma::typeof(1.0u"Å")`: The length scale of the potential.
- `cutoff::Float64`: The cutoff distance for the potential, in units of `sigma`.
- `shift::typeof(0.0u"eV")`: The energy shift applied to the potential, calculated at the cutoff distance.

"""
struct LJParameters <: LennardJonesParametersSets
    epsilon::typeof(1.0u"eV")
    sigma::typeof(1.0u"Å")
    cutoff::Float64
    shift::typeof(0.0u"eV")
end

"""
    LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
A constructor for the LJParameters struct with default values for the 
Lennard-Jones potential with no cutoff or shift. The `shift` parameter can be 
specified as a boolean, if `true`, the shift energy is calculated automatically 
at the cutoff distance; or as a `typeof(0.0u"eV")`, in which case the value is used directly.

# Example

```jldoctest
julia> lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)
LJParameters(0.1 eV, 2.5 Å, 3.5, 0.0 eV)

julia> lj = LJParameters(sigma=2.5)
LJParameters(1.0 eV, 2.5 Å, Inf, 0.0 eV)

julia> lj = LJParameters(cutoff=3.5,shift=5.0)
LJParameters(1.0 eV, 1.0 Å, 3.5, 5.0 eV)

julia> lj = LJParameters(cutoff=3.5,shift=true)
LJParameters(1.0 eV, 1.0 Å, 3.5, -0.0021747803916549904 eV)

julia> lj = LJParameters(cutoff=3.5,shift=false)
LJParameters(1.0 eV, 1.0 Å, 3.5, 0.0 eV)

```
"""
function LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
    if isa(shift,Bool) && shift
        shiftenergy = lj_energy(Float64(epsilon)*u"eV", Float64(sigma)*u"Å", cutoff*Float64(sigma)*u"Å")  
    else
        shiftenergy = Float64(shift)*u"eV"
    end
    return LJParameters(epsilon*u"eV", sigma*u"Å", cutoff, shiftenergy)
end


"""
    lj_energy(epsilon::typeof(1.0u"eV"), sigma::typeof(1.0u"Å"), r::typeof(1.0u"Å"))

Compute the Lennard-Jones potential energy between two particles.

The Lennard-Jones potential energy is given by the equation:

```math
V(r_{ij}) = 4\\varepsilon_{ij} \\left[\\left(\\frac{\\sigma_{ij}}{r_{ij}}\\right)^{12} - \\left(\\frac{\\sigma_{ij}}{r_{ij}}\\right)^{6}\\right]
```

where `epsilon` is the energy scale, `sigma` is the distance scale, and `r` is the distance between the particles.

# Arguments
- `epsilon::typeof(1.0u"eV")`: The energy scale of the potential.
- `sigma::typeof(1.0u"Å")`: The distance scale of the potential.
- `r::typeof(1.0u"Å")`: The distance between the particles.

# Returns
- The Lennard-Jones potential energy between the particles.

"""
function lj_energy(epsilon::typeof(1.0u"eV"), sigma::typeof(1.0u"Å"), r::typeof(1.0u"Å"))
    r6 = (sigma / r)^6
    r12 = r6^2
    return 4 * epsilon * (r12 - r6)
end


"""
    lj_energy(r::typeof(1.0u"Å"), lj::LJParameters)

Compute the Lennard-Jones energy between two particles at a given distance.

# Arguments
- `r::typeof(1.0u"Å")`: The distance between the particles.
- `lj::LJParameters`: The Lennard-Jones parameters.

# Returns
- `0.0u"eV"` if the distance is greater than the cutoff distance.
- The Lennard-Jones energy minus the shift otherwise.
"""
function lj_energy(r::typeof(1.0u"Å"), lj::LJParameters)
    if r > lj.cutoff * lj.sigma
        return 0.0u"eV"
    else
        return lj_energy(lj.epsilon, lj.sigma, r) - lj.shift
    end
end



"""
    struct CompositeLJParameters{C} <: LennardJonesParametersSets

CompositeLJParameters is a struct that represents a set of composite Lennard-Jones parameters.

# Fields
- `lj_param_sets::Matrix{LJParameters}`: A matrix of LJParameters representing the LJ parameter sets.

# Type Parameters
- `C::Int`: The number of composite parameter sets.

"""
struct CompositeLJParameters{C} <: LennardJonesParametersSets
    lj_param_sets::Matrix{LJParameters}
    function CompositeLJParameters{C}(lj_param_sets::Matrix{LJParameters}) where C
        if size(lj_param_sets) != (C, C)
            throw(ArgumentError("the size of the matrix is not compatible with the number of components."))
        end
        new{C}(lj_param_sets)
    end
end

"""
    CompositeLJParameters(c::Int, ljs::Vector{LJParameters})

Construct a `CompositeLJParameters` object from a vector of LJParameters.

# Arguments
- `c::Int`: The number of components.
- `ljs::Vector{LJParameters}`: A vector of LJParameters. 
The number of elements in the vector must be equal to `c^2` or `c*(c+1)/2`. 
The former case is for a full flattened matrix of LJParameters, useful when 
the interactions are asymmetric, i.e., `epsilon_ij != epsilon_ji`. The latter
case is for symmetric interactions, i.e., `epsilon_ij = epsilon_ji`, hence only
the upper triangular part of the matrix is needed.

# Returns
- A `CompositeLJParameters` object.

"""
function CompositeLJParameters(c::Int, ljs::Vector{LJParameters})
    if length(ljs) == c^2
        # If the number of LJParameters sets is equal to the number 
        # of elements in the matrix, assuming that Vector{LJParameters} 
        # is a flattened matrix, then reshape the vector into a matrix.
        return CompositeLJParameters{c}(reshape(ljs, c, c))
    elseif length(ljs) == c*(c+1)/2
        # If the number of LJParameters sets is equal to the number
        # of elements in the upper triangular part of the matrix, then
        # construct the full matrix from the upper triangular part.
        ljmatrix = Matrix{LJParameters}(undef, c, c)
        k = 1
        for i in 1:c
            for j in i:c
            ljmatrix[i, j] = ljs[k]
            if i != j
                ljmatrix[j, i] = ljs[k]
            end
            k += 1
            end
        end
        return CompositeLJParameters{c}(ljmatrix)
    else
        throw(ArgumentError("the number of LJParameters sets is not compatible with the number of components."))
    end
end

end # module Potentials