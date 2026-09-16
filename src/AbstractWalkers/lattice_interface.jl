# ─────────────────────────────────────────────────────────────────────────────
# The AbstractLattice occupancy interface
#
# Sampling code uses these accessors instead of concrete storage fields.
# `MLattice` and `AtomicLattice` implement the same occupation contract.
#
# AtomicLattice supplies neighbor shells in the same site-index space and a
# validated reflection table for periodic point-inversion moves. Reflection is
# deliberately consumed inside `MonteCarloMoves` rather than exposed here: it
# is geometry-specific proposal machinery, not part of the occupancy contract.
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
# AtomicLattice uses the same component-mask representation as MLattice. Its
# constructor and mutation API additionally maintain single-site exclusion,
# because the ASE cache cannot represent two adsorbates at one site.

n_occupied(lattice::AtomicLattice, c::Int=1) = sum(lattice.components[c])
is_occupied(lattice::AtomicLattice, i::Int, c::Int=1) = lattice.components[c][i]
occupied_indices(lattice::AtomicLattice, c::Int=1) = findall(lattice.components[c])
empty_indices(lattice::AtomicLattice, c::Int=1) = findall(.!lattice.components[c])

function set_occupied!(lattice::AtomicLattice, i::Int, v::Bool, c::Int=1)
    checkbounds(lattice.components, c)
    checkbounds(lattice.components[c], i)
    if v
        for component in lattice.components
            component[i] = false
        end
    end
    lattice.components[c][i] = v
    lattice.ase_dirty = true
    return lattice
end

function swap_sites!(lattice::AtomicLattice, a::Int, b::Int)
    a == b && return lattice
    for component in lattice.components
        component[a], component[b] = component[b], component[a]
    end
    lattice.ase_dirty = true
    return lattice
end

function neighbor_shell(lattice::AtomicLattice, i::Int, shell::Int=1)
    return lattice.neighbors[i][shell]
end
