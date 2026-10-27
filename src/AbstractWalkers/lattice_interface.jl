"""
    n_occupied(lattice::AtomicLattice, c::Int=1) -> Int

Number of occupied sites in component `c`.
"""
n_occupied(lattice::AtomicLattice, c::Int=1) = sum(lattice.components[c])

"""
    is_occupied(lattice::AtomicLattice, i::Int, c::Int=1) -> Bool

Whether site `i` of component `c` is occupied.
"""
is_occupied(lattice::AtomicLattice, i::Int, c::Int=1) = lattice.components[c][i]

"""
    occupied_indices(lattice::AtomicLattice, c::Int=1) -> Vector{Int}

Indices of the occupied sites in component `c`.
"""
occupied_indices(lattice::AtomicLattice, c::Int=1) = findall(lattice.components[c])

"""
    empty_indices(lattice::AtomicLattice, c::Int=1) -> Vector{Int}

Indices of the empty sites in component `c`.
"""
empty_indices(lattice::AtomicLattice, c::Int=1) = findall(.!lattice.components[c])

"""
    neighbor_shell(lattice::AtomicLattice, i::Int, shell::Int=1) -> Vector{Int}

Indices of site `i`'s neighbors in the given shell.
"""
neighbor_shell(lattice::AtomicLattice, i::Int, shell::Int=1) = lattice.neighbors[i][shell]

"""
    set_occupied!(lattice::AtomicLattice, i::Int, v::Bool, c::Int=1) -> lattice

Set the occupancy of site `i` in component `c`, returning the lattice.
Setting `v=true` removes any other species at that site, so each site holds
at most one adsorbate.
"""
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

"""
    swap_sites!(lattice::AtomicLattice, a::Int, b::Int) -> lattice

Exchange the occupancy of sites `a` and `b`, across every component.
"""
function swap_sites!(lattice::AtomicLattice, a::Int, b::Int)
    checkbounds(lattice.components[1], a)
    checkbounds(lattice.components[1], b)
    a == b && return lattice
    all(component[a] == component[b] for component in lattice.components) &&
        return lattice
    for component in lattice.components
        component[a], component[b] = component[b], component[a]
    end
    lattice.ase_dirty = true
    return lattice
end

# Sampling masks use direct writes except when an ASE cache must be invalidated.
@inline function _set_lattice_occupied!(lattice::AbstractLattice, i::Int, v::Bool, c::Int=1)
    lattice.components[c][i] = v
    return lattice
end

@inline _set_lattice_occupied!(lattice::AtomicLattice, i::Int, v::Bool, c::Int=1) =
    set_occupied!(lattice, i, v, c)
