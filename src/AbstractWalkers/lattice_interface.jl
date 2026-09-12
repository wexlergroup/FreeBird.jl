# ─────────────────────────────────────────────────────────────────────────────
# The AbstractLattice occupancy interface
#
# Sampling code that works for one lattice type should work for both. Today it
# does not: the grand-canonical sampler, the Monte Carlo moves and the I/O paths
# reach into `lattice.components[1]` directly, which is an `MLattice` field, so
# every one of those call sites is `MLattice`-only by construction rather than
# by intent.
#
# These accessors are the contract that replaces those reads. They are
# deliberately small and deliberately boring: each one is a field access for
# `MLattice` and a field access for `AtomicLattice`, and the point is not that
# they do anything clever but that a caller written against them stops caring
# which type it has.
#
# `neighbor_shell` is not implementable for `AtomicLattice` yet and throws rather
# than guessing. A wrong neighbour list produces a plausible energy, which is the
# worst failure mode available to this package.
#
# The plan's ninth accessor, `reflect_site`, is deliberately absent. Reflection
# through a pivot is arithmetic on a regular grid; its helpers `_site_to_grid`
# and `_reflect_site` live in `MonteCarloMoves`, which loads *after* this module,
# so a method here could not call them. It is also used by exactly one caller,
# the geometric cluster move, which is `MLattice`-only and stays that way until
# `AtomicLattice` can supply a reflection table. Adding it now would buy nothing
# and cost a cross-module dependency in the wrong direction.
# ─────────────────────────────────────────────────────────────────────────────

"""
    n_occupied(lattice::AbstractLattice, c::Int=1) -> Int

Number of occupied sites in component `c`.
"""
function n_occupied end

"""
    is_occupied(lattice::AbstractLattice, i::Int, c::Int=1) -> Bool

Whether site `i` of component `c` is occupied.
"""
function is_occupied end

"""
    set_occupied!(lattice::AbstractLattice, i::Int, v::Bool, c::Int=1) -> lattice

Set the occupancy of site `i` in component `c`, returning the lattice.
"""
function set_occupied! end

"""
    occupied_indices(lattice::AbstractLattice, c::Int=1) -> Vector{Int}

Indices of the occupied sites in component `c`.
"""
function occupied_indices end

"""
    empty_indices(lattice::AbstractLattice, c::Int=1) -> Vector{Int}

Indices of the empty sites in component `c`.
"""
function empty_indices end

"""
    swap_sites!(lattice::AbstractLattice, a::Int, b::Int) -> lattice

Exchange the occupancy of sites `a` and `b`, across every component.
"""
function swap_sites! end

"""
    neighbor_shell(lattice::AbstractLattice, i::Int, shell::Int=1) -> Vector{Int}

Indices of site `i`'s neighbours in the given shell, in the same index space as
the occupancy accessors.

That last clause is the whole content of this function. A neighbour list indexed
over a different set of sites than the occupancy is not a neighbour list, it is
a source of plausible wrong answers.
"""
function neighbor_shell end

# ── MLattice ─────────────────────────────────────────────────────────────────

n_occupied(lattice::MLattice, c::Int=1) = sum(lattice.components[c])
is_occupied(lattice::MLattice, i::Int, c::Int=1) = lattice.components[c][i]
occupied_indices(lattice::MLattice, c::Int=1) = findall(lattice.components[c])
empty_indices(lattice::MLattice, c::Int=1) = findall(.!lattice.components[c])

function set_occupied!(lattice::MLattice, i::Int, v::Bool, c::Int=1)
    lattice.components[c][i] = v
    return lattice
end

function swap_sites!(lattice::MLattice{C,G}, a::Int, b::Int) where {C,G}
    a == b && return lattice
    for comp in 1:C
        lattice.components[comp][a], lattice.components[comp][b] =
            lattice.components[comp][b], lattice.components[comp][a]
    end
    return lattice
end

neighbor_shell(lattice::MLattice, i::Int, shell::Int=1) = lattice.neighbors[i][shell]

# ── AtomicLattice ────────────────────────────────────────────────────────────
#
# Occupancy lives in a single `Vector{Bool}` over `all_sites`, so the component
# argument exists only to match the interface and must be 1.

@inline function _check_single_component(c::Int)
    c == 1 || throw(ArgumentError(
        "AtomicLattice carries one occupancy mask over all_sites; got component $c"))
    return nothing
end

function n_occupied(lattice::AtomicLattice, c::Int=1)
    _check_single_component(c)
    return sum(lattice.occupations)
end

function is_occupied(lattice::AtomicLattice, i::Int, c::Int=1)
    _check_single_component(c)
    return lattice.occupations[i]
end

function occupied_indices(lattice::AtomicLattice, c::Int=1)
    _check_single_component(c)
    return findall(lattice.occupations)
end

function empty_indices(lattice::AtomicLattice, c::Int=1)
    _check_single_component(c)
    return findall(.!lattice.occupations)
end

function set_occupied!(lattice::AtomicLattice, i::Int, v::Bool, c::Int=1)
    _check_single_component(c)
    lattice.occupations[i] = v
    lattice.ase_dirty = true
    return lattice
end

function swap_sites!(lattice::AtomicLattice, a::Int, b::Int)
    a == b && return lattice
    lattice.occupations[a], lattice.occupations[b] =
        lattice.occupations[b], lattice.occupations[a]
    lattice.ase_dirty = true
    return lattice
end

function neighbor_shell(lattice::AtomicLattice, i::Int, shell::Int=1)
    throw(ArgumentError(
        "neighbor_shell is not available for AtomicLattice yet. Its `neighbors` " *
        "field is indexed over the substrate grid built by get_lattice_positions, " *
        "while occupancy is indexed over `all_sites` — the union of ontop, bridge " *
        "and hollow positions selected by type_of_sites. The two are different " *
        "site sets in a different order, so returning neighbors[i] here would " *
        "hand back a neighbour list for a site other than the one asked about. " *
        "Re-indexing the neighbour list onto all_sites is the outstanding work."))
end
