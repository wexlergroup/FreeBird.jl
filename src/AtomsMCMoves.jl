module AtomsMCMoves

using ExtXYZ
using AtomsBase
using Setfield
using Distributions
using Unitful
using StaticArrays

using ..Potentials
using ..AbstractWalkers
using ..EnergyEval

export periodic_boundary_wrap!
export MC_random_walk!, MC_nve_walk!


function periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T
    pbc = system.boundary_conditions
    box = system.bounding_box
    pos = SVector{3,T}([pbc[i]==Periodic() ? mod(pos[i], box[i][i]) : pos[i] for i in eachindex(pbc)])
    return pos
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
        deleteat!(config.position, i_at)
        insert!(config.position, i_at, pos)
        energy = interaction_energy(config, lj; frozen=frozen) + e_shift
        if energy >= emax
            # reject the move, revert to original position
            deleteat!(config.position, i_at)
            insert!(config.position, i_at, orig_pos)
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
                        demon_energy_threshold=Inf*u"eV"
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
    deleteat!(config.position, i_at)
    insert!(config.position, i_at, pos)
    new_energy = interaction_energy(config, lj; frozen=frozen) + e_shift
    ΔE = new_energy - orig_energy
    if ΔE <= 0.0u"eV" && e_demon < demon_energy_threshold
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
        deleteat!(config.position, i_at)
        insert!(config.position, i_at, orig_pos)
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon has no enough energy, reject"
    end
    return accept, at, e_demon
end

function additional_demon_walk!(
                    e_demon::typeof(0.0u"eV"),
                    at::AtomWalker,
                    lj::LJParameters, 
                    step_size::Float64;
                    e_demon_tolerance=1e-10u"eV",
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
            @info "Further reduced step size to $(step_size/10)."
            allow_step_size_reduction = false
        end
    end
    if additional_steps >= max_add_steps
        @warn "Maximum additional steps reached, breaking. NVE condition not satisfied with an energy by a factor of $(e_demon/e_demon_tolerance)."
    else
        accept_this_walker = true
    end
    @info "Performed $additional_steps additional demon walk.
    New accept rate: $(new_accept_count/additional_steps).
    Final demon energy: $e_demon."
    return accept_this_walker, new_accept_count/additional_steps, at
end

function MC_nve_walk!(
                    n_steps::Int, 
                    at::AtomWalker, 
                    lj::LJParameters, 
                    step_size::Float64; 
                    e_demon_tolerance=1e-10u"eV",
                    demon_energy_threshold=Inf*u"eV",
                    max_add_steps::Int=1_000_000)
    accept_this_walker = false
    e_demon = 0.0u"eV"
    accept_count::Int64 = 0
    initial_energy  = at.energy
    @info "Initial energy: $(at.energy)"
    demon_energies = Array{typeof(0.0u"eV")}(undef, n_steps)
    for i in 1:n_steps
        accept, at, e_demon = single_atom_demon_walk!(at, lj, step_size; 
                                e_demon=e_demon, 
                                demon_energy_threshold=demon_energy_threshold,
                                )
        demon_energies[i] = e_demon
        if accept
            accept_count += 1
        end
    end
    #TODO: statistics on demon energies
    if e_demon > e_demon_tolerance
        @info "Demon has too much energy remain (>tolerance=$(e_demon_tolerance)): $e_demon, performing more demon walk.
        Current accept rate: $(accept_count/n_steps). Step size is lowered."
        step_size /= 10
        accept_this_walker, _, at = additional_demon_walk!(e_demon, at, lj, step_size; 
                                                e_demon_tolerance=e_demon_tolerance, 
                                                max_add_steps=max_add_steps)
    else
        accept_this_walker = true
    end

    final_energy = at.energy + e_demon
    @info "Final energy: $(at.energy)"
    @debug "initial_energy: ", initial_energy, " final_energy: ", final_energy
    return accept_this_walker, accept_count/n_steps, at
end


end # module AtomsMonteCarloMoves