#=
This file contains energy evaluation functions for atomistic systems using many-body potentials.
Currently, it includes implementations for the Gupta potential. In principle, it should also work 
for other many-body potentials that can be defined using the `two_body_energy` and `many_body_energy` functions.
=#

"""
    intra_component_energy(at::AbstractSystem, pot::GuptaParameters)

Calculate the energy from interactions between particles using the Gupta potential for a single component.
"""
function intra_component_energy(at::AbstractSystem, pot::GuptaParameters)
    # Calculate the energy from interactions between particles using the Gupta potential.
    # The energy is calculated by summing the pairwise interactions between the free particles.

    pair_energies = zeros(Float64, length(at))
    bond_energies = zeros(Float64, length(at))

    Threads.@threads for i in 1:length(at)
        for j in 1:length(at)
            if i != j
                r = pbc_dist(position(at, i), position(at, j), at)
                bond_energies[i] += ustrip(many_body_energy(r, pot))
                pair_energies[i] += ustrip(two_body_energy(r, pot))
            end
        end
    end
    return total_energy(pair_energies, bond_energies, pot)
end

"""
    interacting_energy(at::AbstractSystem, 
                            pot::CompositeParameterSets{C,GuptaParameters},
                            list_num_par::Vector{Int}
                            ) where C
Calculate the energy from interactions between particles using the Gupta potential for a multi-component system.
"""
function interacting_energy(at::AbstractSystem, 
                            pot::CompositeParameterSets{C,GuptaParameters},
                            list_num_par::Vector{Int}
                            ) where C
    # Calculate the energy from interactions between particles using the Gupta potential.
    # The energy is calculated by summing the pairwise interactions between the free particles.

    pair_energies = zeros(Float64, length(at))
    bond_energies = zeros(Float64, length(at))

    # Expand list_num_par so that each element n is repeated n times
    expanded_list_num_par::Vector{Int} = vcat([fill(ind, x) for (ind, x) in enumerate(list_num_par)]...)
    # println("expanded_list_num_par: $expanded_list_num_par") # DEBUG

    Threads.@threads for i in 1:length(at)
        for j in 1:length(at)
            if i != j
                element_a, element_b = expanded_list_num_par[i], expanded_list_num_par[j]
                # println("element_a: $element_a, element_b: $element_b") # DEBUG
                r = pbc_dist(position(at, i), position(at, j), at)
                bond_energies[i] += ustrip(many_body_energy(r, pot.param_sets[element_a, element_b]))
                pair_energies[i] += ustrip(two_body_energy(r, pot.param_sets[element_a, element_b]))
            end
        end
    end
    return total_energy(pair_energies, bond_energies, pot)
end


"""
    interacting_energy(system::AbstractSystem, calc::ASECalculator)

Evaluate total potential energy (eV) of `system` using the wrapped ASE
calculator.
"""
function interacting_energy(system::AtomsBase.AbstractSystem, calc::PyCalculator)
    return AtomsCalculators.potential_energy(system, calc.calc)
end

# variant with list_num_par/frozen passed by callers in LJ workflow
function interacting_energy(system::AtomsBase.AbstractSystem,
                            calc::PyCalculator,
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool})
    return interacting_energy(system, calc)
end