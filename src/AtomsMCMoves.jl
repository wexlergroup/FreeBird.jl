module AtomsMCMoves

using ExtXYZ
using AtomsBase
using Setfield
using Distributions
using Unitful
using StaticArrays
# using LinearAlgebra

using ..Potentials
using ..AbstractWalkers
using ..EnergyEval

export periodic_boundary_wrap!
export MC_random_walk!, MC_nve_walk!


function periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T
    pbc = system.boundary_conditions
    box = system.bounding_box
    new_pos = Vector{typeof(0.0u"Å")}(undef, 3)
    for i in eachindex(pbc)
        if pbc[i] == Periodic()
            new_pos[i] = mod(pos[i], box[i][i])
        else
            new_pos[i] = pos[i]
        end
    end
    pos = SVector{3,T}(new_pos)
    return pos
end

function mean_sq_displacement(at::AtomWalker, at_orig::AtomWalker)
    distsq::typeof(0.0u"Å"^2) = 0.0u"Å"^2
    frozen = at.num_frozen_part
    for i in (frozen+1):length(at.configuration)
        dist::typeof(0.0u"Å") = pbc_dist(at.configuration.position[i], at_orig.configuration.position[i], at.configuration)
        distsq += dist^2
    end
    return distsq/(length(at.configuration)-frozen)
end


function single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64) where T
    (dx, dy, dz) = (rand(Uniform(-step_size,step_size)) for _ in 1:3) .* unit(T)
    pos = pos .+ (dx, dy, dz)
    return pos 
end


"""
    random_walk(n_steps::Int, at::Atoms, lj::LJParameters, step_size::Float64, emax::Float64, frozen::Int)
Perform a random walk of `n_steps` steps on the atoms in `at` using a step size of `step_size`.
"""
function MC_random_walk!(
                    n_steps::Int, 
                    at::AtomWalker, 
                    lj::LJParameters, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV")
                    )
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        frozen = at.num_frozen_part
        e_shift = at.energy_frozen_part
        i_at = rand((frozen+1):length(config))
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        energy = interaction_energy(config, lj; frozen=frozen) + e_shift
        if energy >= emax
            # reject the move, revert to original position
            config.position[i_at] = orig_pos
        else
            at.energy = energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end




function single_atom_demon_walk!(
                        at::AtomWalker, 
                        lj::LJParameters, 
                        step_size::Float64;
                        e_demon=0.0u"eV",
                        demon_energy_threshold=Inf*u"eV",
                        demon_gain_threshold=Inf*u"eV",
                        )
    accept = false
    orig_energy = at.energy
    config = at.configuration
    frozen = at.num_frozen_part
    e_shift = at.energy_frozen_part
    i_at = rand((frozen+1):length(config))
    pos = position(config, i_at)
    orig_pos = deepcopy(pos)
    pos = single_atom_random_walk!(pos, step_size)
    pos = periodic_boundary_wrap!(pos, config)
    config.position[i_at] = pos
    new_energy = interaction_energy(config, lj; frozen=frozen) + e_shift
    ΔE = new_energy - orig_energy
    if -demon_gain_threshold <= ΔE <= 0.0u"eV" && e_demon-ΔE < demon_energy_threshold
        e_demon -= ΔE
        accept = true
        at.energy = new_energy
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon gains energy, accept"
    elseif 0.0u"eV" < ΔE <= e_demon
        e_demon -= ΔE
        accept = true
        at.energy = new_energy
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon gives energy, accept"
    else
        config.position[i_at] = orig_pos
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon has no enough energy, reject"
    end
    return accept, at, e_demon
end

function additional_demon_walk!(
                    e_demon::typeof(0.0u"eV"),
                    at::AtomWalker,
                    lj::LJParameters, 
                    step_size::Float64;
                    e_demon_tolerance=1e-9u"eV",
                    max_add_steps::Int=1_000_000
                    )
    accept_this_walker = false
    new_accept_count::Int = 0
    additional_steps::Int = 1
    allow_step_size_reduction = true
    while e_demon > e_demon_tolerance && additional_steps < max_add_steps
        accept, at, e_demon = single_atom_demon_walk!(at, lj, step_size; e_demon=e_demon, demon_energy_threshold=0.0u"eV")
        if accept
            new_accept_count += 1
        end
        additional_steps += 1
        if additional_steps > (max_add_steps ÷ 2) && allow_step_size_reduction
            step_size /= 10
            @info "Further reduced step size to $(step_size)."
            allow_step_size_reduction = false
        end
    end
    if additional_steps >= max_add_steps
        @warn "Maximum additional steps reached NVE condition not satisfied with an energy by a factor of $(e_demon/e_demon_tolerance)."
    else
        accept_this_walker = true
    end
    @info "Performed $additional_steps additional demon walk."
    @debug "New accept rate: $(new_accept_count/additional_steps)."
    @debug "Final demon energy: $e_demon."
    return accept_this_walker, new_accept_count/additional_steps, at, e_demon
end

function MC_nve_walk!(
                    n_steps::Int, 
                    at::AtomWalker, 
                    lj::LJParameters, 
                    step_size::Float64; 
                    e_demon_tolerance=1e-9u"eV",
                    demon_energy_threshold=Inf*u"eV",
                    demon_gain_threshold=Inf*u"eV",
                    max_add_steps::Int=1_000_000)
    accept_this_walker = false
    e_demon = 0.0u"eV"
    at_original = deepcopy(at)
    accept_count::Int64 = 0
    initial_energy  = at.energy
    @debug "Initial energy: $(at.energy)"
    demon_energies = Array{typeof(0.0u"eV")}(undef, n_steps)
    for i in 1:n_steps
        accept, at, e_demon = single_atom_demon_walk!(at, lj, step_size; 
                                e_demon=e_demon, 
                                demon_energy_threshold=demon_energy_threshold,
                                demon_gain_threshold=demon_gain_threshold,
                                )
        demon_energies[i] = e_demon
        if accept
            accept_count += 1
        end
    end
    temp_estimate = mean(demon_energies)
    temp_per_K = temp_estimate/(8.61733e-5u"eV/K")
    @info "Estimated temperature: $temp_estimate, or $temp_per_K."
    if accept_count == 0
        @warn "No demon walk was accepted, reverting to original configuration."
    elseif e_demon > e_demon_tolerance
        step_size /= 10
        @info "Demon has too much energy remain (>=$e_demon_tolerance): $e_demon, performing more demon walk. Step size is lowered to $step_size."
        accept_this_walker, _, at, e_demon = additional_demon_walk!(e_demon, at, lj, step_size; 
                                                e_demon_tolerance=e_demon_tolerance, 
                                                max_add_steps=max_add_steps)
    else
        accept_this_walker = true
    end

    final_energy = at.energy + e_demon
    @debug "Final energy: $(at.energy)"
    @debug "initial_energy: $initial_energy, final_energy: $final_energy"
    if !(final_energy ≈ initial_energy)
        @error "Energy conservation not satisfied: initial energy: $initial_energy, final energy: $final_energy."
    end

    if !accept_this_walker
        return accept_this_walker, accept_count/n_steps, at_original, demon_energies, temp_estimate
    else
        mean_sq_displace = mean_sq_displacement(at, at_original)
        @info "Mean squared displacement: $mean_sq_displace."
        return accept_this_walker, accept_count/n_steps, at, demon_energies, temp_estimate
    end
end


end # module AtomsMonteCarloMoves