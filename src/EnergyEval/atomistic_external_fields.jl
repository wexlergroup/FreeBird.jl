# Energy evaluators for `ExternalFieldPotential`: the wrapped pair potential's value
# plus the one-body field term. No existing method body or signature changes; these
# methods are more specific than the `SingleComponentPotential{Pairwise}` ones.

# Sum of the field energy over the atoms of the components selected by `want_frozen`
# (the free components when false, the frozen ones when true).
function _field_energy_sum(field, at::AbstractSystem, list_num_par::Vector{Int},
                           frozen::Vector{Bool}, want_frozen::Bool)
    energy = 0.0u"eV"
    offset = 0
    for (c, n) in enumerate(list_num_par)
        if frozen[c] == want_frozen
            for i in (offset + 1):(offset + n)
                energy += external_energy(field, position(at, i), at)
            end
        end
        offset += n
    end
    return energy
end

"""
    single_site_energy(index::Int, at::AbstractSystem, pot::ExternalFieldPotential, list_num_par::Vector{Int})

The energy of one atom under an `ExternalFieldPotential`: the wrapped pair potential's
single-site energy plus the field energy at the atom's position. The field is evaluated
first, and `+Inf` is returned before the O(N) pair sum when the atom lies outside the
accessible region; no random draw is involved, so the short-circuit is stream-neutral.
"""
function single_site_energy(index::Int,
                            at::AbstractSystem,
                            pot::ExternalFieldPotential,
                            list_num_par::Vector{Int})
    u_ext = external_energy(pot.field, position(at, index), at)
    isfinite(u_ext) || return u_ext
    return single_site_energy(index, at, pot.pair, list_num_par) + u_ext
end

"""
    interacting_energy(at::AbstractSystem, pot::ExternalFieldPotential, list_num_par::Vector{Int}, frozen::Vector{Bool})

The interacting energy under an `ExternalFieldPotential`: the wrapped pair potential's
free-free and free-frozen energy plus the field energy of the free atoms. The field
energy of the frozen atoms is a constant carried by `frozen_energy`.
"""
function interacting_energy(at::AbstractSystem,
                            pot::ExternalFieldPotential,
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool})
    e_pair = interacting_energy(at, pot.pair, list_num_par, frozen)
    return e_pair + _field_energy_sum(pot.field, at, list_num_par, frozen, false)
end

"""
    interacting_energy(at::AbstractSystem, pot::ExternalFieldPotential)

The total energy under an `ExternalFieldPotential` with every atom treated as free: the
wrapped pair potential's energy plus the field energy of every atom.
"""
function interacting_energy(at::AbstractSystem, pot::ExternalFieldPotential)
    e_ext = 0.0u"eV"
    for i in 1:length(at)
        e_ext += external_energy(pot.field, position(at, i), at)
    end
    return interacting_energy(at, pot.pair) + e_ext
end

"""
    frozen_energy(at::AbstractSystem, pot::ExternalFieldPotential, list_num_par::Vector{Int}, frozen::Vector{Bool})

The frozen energy under an `ExternalFieldPotential`: the wrapped pair potential's
frozen-frozen energy plus the field energy of the frozen atoms.
"""
function frozen_energy(at::AbstractSystem,
                       pot::ExternalFieldPotential,
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool})
    e_pair = frozen_energy(at, pot.pair, list_num_par, frozen)
    return e_pair + _field_energy_sum(pot.field, at, list_num_par, frozen, true)
end
