"""
    inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::AbstractPotential)

Compute the energy between two components of a system using a specified (pairwise) potential.

# Arguments
- `at1::AbstractSystem`: The first component of the system.
- `at2::AbstractSystem`: The second component of the system.
- `pot::AbstractPotential`: The potential used to compute the energy.

# Returns
- `energy`: The energy between the two components.

"""
function inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::SingleComponentPotential{Pairwise})
    # build pairs of particles
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2)]
    # @show pairs # DEBUG
    energy = Array{typeof(0.0u"eV"), 1}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        # @show i,j # DEBUG
        (i, j) = pairs[k]
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energy[k] = pair_energy(r, pot)
    end
    # energy = energy*u"eV"
    return sum(energy)
end

"""
    intra_component_energy(at::AbstractSystem, pot::AbstractPotential)

Compute the energy within a component of a system using a specified (pairwise) potential.

# Arguments
- `at::AbstractSystem`: The component of the system.
- `pot::AbstractPotential`: The potential used to compute the energy.

# Returns
- `energy`: The energy within the component.

"""
function intra_component_energy(at::AbstractSystem, pot::SingleComponentPotential{Pairwise}) 
    # num_pairs = length(at) * (length(at) - 1) ÷ 2
    pairs = Array{Tuple{Int,Int}, 1}()
    for i in 1:length(at)
        for j in (i+1):length(at)
            push!(pairs, (i, j))
        end
    end
    # @info "num_pairs: $num_pairs, length(pairs): $(length(pairs))"
    energies = Vector{typeof(0.0u"eV")}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        (i, j) = pairs[k]
        r = pbc_dist(position(at, i), position(at, j), at)
        energies[k] = pair_energy(r, pot)
        # @info "interacting pair: [$(i),$(j)] $(lj_energy(r,lj))"
    end
    return sum(energies)
end
