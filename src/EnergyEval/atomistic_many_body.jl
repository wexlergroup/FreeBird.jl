function interacting_energy(at::AbstractSystem, pot::SingleComponentPotential{ManyBody})
    # Calculate the energy from interactions between particles using the a single-component potential.
    # The energy is calculated by summing the pairwise interactions between the free particles.
    
    # Build pairs of particles without double counting
    pairs = Array{Tuple{Int,Int}, 1}()
    for i in 1:length(at)
        for j in (i+1):length(at)
            push!(pairs, (i, j))
        end
    end
    
    energies_two_body = zeros(Float64, length(pairs))
    energies_many_body = zeros(Float64, length(at))
    
    Threads.@threads for k in eachindex(pairs)
        (i, j) = pairs[k]
        r = pbc_dist(position(at, i), position(at, j), at)
        energies_two_body[k] = ustrip(two_body_energy(r, pot))
    end
    
    Threads.@threads for i in eachindex(at.position)
        for j in eachindex(at.position)
            if i != j
                r = pbc_dist(position(at, i), position(at, j), at)
                energies_many_body[i] += ustrip(many_body_energy(r, pot))
            end
        end
    end
    # @show energies_two_body, energies_many_body # Debugging output to check energies
    return total_energy(energies_two_body, energies_many_body, pot)
end