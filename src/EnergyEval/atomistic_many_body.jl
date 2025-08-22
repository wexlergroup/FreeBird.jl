function interacting_energy(at::AbstractSystem, pot::SingleComponentPotential{ManyBody})
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