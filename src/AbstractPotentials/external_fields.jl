# One-body external fields and the wrapper potential that adds them to a
# single-component pairwise potential.

"""
    AbstractExternalField

Abstract type for one-body external fields: an energy that depends only on each atom's
position, such as a smooth wall confining a fluid in a pore or a uniform field across a
slit. A concrete field must be serializable (no closures: the distributed canonical steps
send the potential to worker processes) and implement three methods:

- `external_energy(field, pos, cell)`: the field energy of one atom at `pos`, a
  `typeof(0.0u"eV")` that is `+Inf` outside the accessible region and never `NaN`;
- `accessible(field, pos, cell)`: whether `pos` lies in the accessible region;
- `accessible_volume(field, cell)`: the volume of the accessible region inside the cell.

`cell` is the simulation system (any `AtomsBase.AbstractSystem`), from which the cell
vectors are read. Fields assume an orthorhombic cell, like `pbc_dist`.
"""
abstract type AbstractExternalField end

"""
    external_energy(field::AbstractExternalField, pos, cell::AbstractSystem)

The one-body energy of an atom at position `pos` (three lengths) in the field: a
`typeof(0.0u"eV")` that is `+Inf` outside the field's accessible region and never `NaN`.
"""
function external_energy end

"""
    accessible(field::AbstractExternalField, pos, cell::AbstractSystem)

Whether the position `pos` lies in the field's accessible region, where its energy is
finite. The restricted initializers draw positions by rejection against it.
"""
function accessible end

"""
    accessible_volume(field::AbstractExternalField, cell::AbstractSystem)

The volume of the field's accessible region inside the orthorhombic cell of `cell`, a
`typeof(1.0u"Å^3")`. Pass it as `V` to the grand-canonical reductions of walkers drawn by
the restricted initializers.
"""
function accessible_volume end

function _orthorhombic_box(cell::AbstractSystem)
    cellv = cell_vectors(cell)
    for i in 1:3, j in 1:3
        if i != j && !iszero(ustrip(cellv[i][j]))
            throw(ArgumentError("external fields assume an orthorhombic cell; found a nonzero off-diagonal cell component"))
        end
    end
    return (cellv[1][1], cellv[2][2], cellv[3][3])
end

function _check_axis(axis::Int)
    1 <= axis <= 3 || throw(ArgumentError("axis must be 1, 2 or 3; got $axis"))
    return axis
end

function _check_table(x::AbstractVector, U::AbstractVector, name::String)
    length(x) == length(U) || throw(ArgumentError("the $name nodes and the energies must have equal length"))
    length(x) >= 2 || throw(ArgumentError("a tabulated field needs at least two nodes"))
    all(diff(ustrip.(u"Å", x)) .> 0.0) || throw(ArgumentError("the $name nodes must be strictly increasing"))
    all(isfinite, ustrip.(u"eV", U)) || throw(ArgumentError("the tabulated energies must be finite; the region outside the table is already +Inf"))
    return nothing
end

# Piecewise-linear interpolation on strictly increasing nodes; the caller has
# already checked x ∈ [xs[1], xs[end]].
function _interpolate(xs::Vector{typeof(1.0u"Å")}, Us::Vector{typeof(0.0u"eV")}, x)
    k = searchsortedlast(xs, x)
    k >= length(xs) && return Us[end]
    t = (x - xs[k]) / (xs[k+1] - xs[k])
    return Us[k] + t * (Us[k+1] - Us[k])
end

"""
    ZeroField()

The field that is zero everywhere; its accessible region is the whole cell. Wrapping a
potential with it reproduces the unwrapped potential's energies exactly.
"""
struct ZeroField <: AbstractExternalField end

external_energy(::ZeroField, pos, cell::AbstractSystem) = 0.0u"eV"
accessible(::ZeroField, pos, cell::AbstractSystem) = true
function accessible_volume(::ZeroField, cell::AbstractSystem)
    box = _orthorhombic_box(cell)
    return box[1] * box[2] * box[3]
end

"""
    TabulatedPlanarField(axis, z, U)

A field that depends on one Cartesian coordinate: `U` tabulated on the strictly
increasing nodes `z` along `axis` (1, 2 or 3), interpolated piecewise-linearly, and
`+Inf` outside `[z[1], z[end]]`. The accessible volume is the cell's cross-section
normal to `axis` times the table's span clipped to the cell.

# Arguments
- `axis::Int`: the Cartesian axis the field varies along.
- `z::AbstractVector`: the nodes, lengths (converted to Å).
- `U::AbstractVector`: the energies at the nodes (converted to eV), all finite.
"""
struct TabulatedPlanarField <: AbstractExternalField
    axis::Int
    z::Vector{typeof(1.0u"Å")}
    U::Vector{typeof(0.0u"eV")}
    function TabulatedPlanarField(axis::Int, z::AbstractVector, U::AbstractVector)
        _check_axis(axis)
        _check_table(z, U, "planar")
        return new(axis, Vector{typeof(1.0u"Å")}(uconvert.(u"Å", z)),
                   Vector{typeof(0.0u"eV")}(uconvert.(u"eV", U)))
    end
end

function external_energy(f::TabulatedPlanarField, pos, cell::AbstractSystem)
    x = pos[f.axis]
    (x < f.z[1] || x > f.z[end]) && return Inf * u"eV"
    return _interpolate(f.z, f.U, x)
end

function accessible(f::TabulatedPlanarField, pos, cell::AbstractSystem)
    x = pos[f.axis]
    return f.z[1] <= x <= f.z[end]
end

function accessible_volume(f::TabulatedPlanarField, cell::AbstractSystem)
    box = _orthorhombic_box(cell)
    lo = max(f.z[1], zero(f.z[1]))
    hi = min(f.z[end], box[f.axis])
    span = max(hi - lo, zero(lo))
    # the product in the cell's own order, so a table spanning the cell gives V_cell bit for bit
    return prod(i == f.axis ? span : box[i] for i in 1:3)
end

"""
    TabulatedRadialField(axis, center, r, U)

A field that depends on the distance ρ from the line parallel to `axis` through
`center`: `U` tabulated on the strictly increasing nodes `r` (with `r[1] = 0`),
interpolated piecewise-linearly, and `+Inf` for ρ > `r[end]`. The distance is measured
directly, without minimum image, so the disc of radius `r[end]` must fit inside the
cell's cross-section; `accessible_volume` (called by the initializers) throws an
`ArgumentError` when it does not. The accessible volume is π `r[end]`² times the cell
length along `axis`. A hard cylinder is a table of zeros.

# Arguments
- `axis::Int`: the Cartesian axis the line runs along.
- `center::AbstractVector`: a point on the line, three lengths (converted to Å); its
  component along `axis` is ignored.
- `r::AbstractVector`: the radial nodes, lengths (converted to Å), starting at zero.
- `U::AbstractVector`: the energies at the nodes (converted to eV), all finite.
"""
struct TabulatedRadialField <: AbstractExternalField
    axis::Int
    center::SVector{3, typeof(1.0u"Å")}
    r::Vector{typeof(1.0u"Å")}
    U::Vector{typeof(0.0u"eV")}
    function TabulatedRadialField(axis::Int, center::AbstractVector, r::AbstractVector, U::AbstractVector)
        _check_axis(axis)
        length(center) == 3 || throw(ArgumentError("center must have three components"))
        _check_table(r, U, "radial")
        iszero(ustrip(u"Å", r[1])) || throw(ArgumentError("the radial nodes must start at r = 0"))
        return new(axis, SVector{3, typeof(1.0u"Å")}(uconvert.(u"Å", center)),
                   Vector{typeof(1.0u"Å")}(uconvert.(u"Å", r)),
                   Vector{typeof(0.0u"eV")}(uconvert.(u"eV", U)))
    end
end

function _radial_distance(f::TabulatedRadialField, pos)
    s = zero(f.r[1])^2
    for i in 1:3
        i == f.axis && continue
        s += (pos[i] - f.center[i])^2
    end
    return sqrt(s)
end

function external_energy(f::TabulatedRadialField, pos, cell::AbstractSystem)
    rho = _radial_distance(f, pos)
    rho > f.r[end] && return Inf * u"eV"
    return _interpolate(f.r, f.U, rho)
end

accessible(f::TabulatedRadialField, pos, cell::AbstractSystem) = _radial_distance(f, pos) <= f.r[end]

function accessible_volume(f::TabulatedRadialField, cell::AbstractSystem)
    box = _orthorhombic_box(cell)
    for i in 1:3
        i == f.axis && continue
        if f.center[i] - f.r[end] < zero(f.r[end]) || f.center[i] + f.r[end] > box[i]
            throw(ArgumentError("the disc of radius r[end] = $(f.r[end]) around the field's axis does not fit the cell's cross-section along axis $i"))
        end
    end
    return pi * f.r[end]^2 * box[f.axis]
end

"""
    ExternalFieldPotential(pair, field)

A single-component pairwise potential plus a one-body external field: the energy of a
configuration is the pair energy of `pair` plus the sum of `external_energy(field, ...)`
over the atoms. The wrapper reaches the grand-canonical nested-sampling kernel, the µVT
kernel and driver, and the canonical walk through its own methods of `single_site_energy`,
`interacting_energy`, `frozen_energy` and `MC_random_walk!`; `pair_energy` is deliberately
not defined for it, so a code path that bypasses those evaluators fails with a
`MethodError` instead of dropping the field silently.

Hard regions: a field that is `+Inf` outside an accessible region targets the reference
measure restricted to that region. The kernels are unchanged (a proposal outside the region
carries `+Inf` and fails), the grand-canonical initializer draws N ~ Poisson(z₀V_acc) with
positions uniform in the region, and the reduction is called with
`V = accessible_volume(field, cell)`. No walker ever holds `+Inf`. The µVT driver needs no
restricted initializer: its activity volume stays `z V_cell` (its kernel inserts uniformly in
the cell, and insertions outside the region carry a zero Boltzmann factor).

# Fields
- `pair::SingleComponentPotential{Pairwise}`: the pair potential (not itself a wrapper).
- `field::AbstractExternalField`: the one-body field.
"""
struct ExternalFieldPotential{P<:SingleComponentPotential{Pairwise}, F<:AbstractExternalField} <: SingleComponentPotential{Pairwise}
    pair::P
    field::F
    function ExternalFieldPotential(pair::P, field::F) where {P<:SingleComponentPotential{Pairwise}, F<:AbstractExternalField}
        pair isa ExternalFieldPotential && throw(ArgumentError("the pair potential of an ExternalFieldPotential cannot itself be a wrapper; combine the fields into one field"))
        return new{P,F}(pair, field)
    end
end

_max_interaction_range(w::ExternalFieldPotential) = _max_interaction_range(w.pair)
