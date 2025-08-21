function inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::GuptaParameters)
    # build pairs of particles
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2)]
    # @show pairs # DEBUG
    energies_rep = Array{typeof(0.0u"eV"), 1}(undef, length(pairs))
    energies_att = Array{typeof(0.0u"eV^2"), 1}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        # @show i,j # DEBUG
        (i, j) = pairs[k]
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energies_rep[k] = gupta_repulsion(r,pot)
        energies_att[k] = gupta_attraction_squared(r,pot)
    end
    total_rep = sum(energies_rep)
    total_att = sum(energies_att)
    total_energy = total_rep - sqrt(total_att) # Gupta potential energy
    return total_energy
end


function intra_component_energy(at::AbstractSystem, pot::GuptaParameters)
    # num_pairs = length(at) * (length(at) - 1) ÷ 2
    pairs = Array{Tuple{Int,Int}, 1}()
    for i in 1:length(at)
        for j in (i+1):length(at)
            push!(pairs, (i, j))
        end
    end
    # @info "num_pairs: $num_pairs, length(pairs): $(length(pairs))"
    energies_rep = Vector{typeof(0.0u"eV")}(undef, length(pairs))
    energies_att = Vector{typeof(0.0u"eV^2")}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        (i, j) = pairs[k]
        r = pbc_dist(position(at, i), position(at, j), at)
        energies_rep[k] = gupta_repulsion(r,pot)
        energies_att[k] = gupta_attraction_squared(r,pot)
        # @info "interacting pair: [$(i),$(j)] $(lj_energy(r,lj))"
    end
    total_rep = sum(energies_rep)
    total_att = sum(energies_att)
    total_energy = total_rep - sqrt(total_att) # Gupta potential energy
    return total_energy
end