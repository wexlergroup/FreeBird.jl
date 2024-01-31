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


function periodic_boundary_wrap!(pos::Vector{T}, system::AbstractSystem) where T
    pbc = system.boundary_conditions
    box = system.bounding_box
    for i in eachindex(pos)
        if pbc[i] == Periodic()
            pos[i] = mod(pos[i], box[i][i])
        end 
    end
    return pos
end



function single_atom_random_walk!(pos::Vector{T}, step_size::Float64) where T
    (dx, dy, dz) = (rand(Uniform(-step_size,step_size)) for _ in 1:3) .* unit(T)
    pos .+= (dx, dy, dz)
    return pos 
end


"""
    random_walk(n_steps::Int, at::Atoms, lj::LJParameters, step_size::Float64, emax::Float64, frozen::Int)
Perform a random walk of `n_steps` steps on the atoms in `at` using a step size of `step_size`.
"""
function MC_random_walk!(
                    n_steps::Int, 
                    at::AbstractSystem, 
                    lj::LJParameters, 
                    step_size::Float64, 
                    emax::Float64; 
                    frozen::Int=0, 
                    e_shift::Float64=0.0
                    )
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        i_at = rand((frozen+1):length(at))
        pos = position(at, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        energy = interaction_energy(at, lj; frozen=frozen) + e_shift
        if energy >= emax
            # reject the move, revert to original position
            at.particles[i_at][:position] .= orig_pos
        else
            at.data[:energy] = energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept, at
end




function single_atom_demon_walk!(
                        at::Atoms, 
                        lj::LJParameters, 
                        step_size::Float64; 
                        frozen::Int=0, 
                        e_shift::Float64=0.0, 
                        e_demon::Float64=0.0,
                        demon_energy_threshold::Float64=0.0
                        )
    accept = false
    i_at = rand((frozen+1):length(at))
    orig_energy = at.system_data.energy
    pos = at.atom_data.position[i_at]
    orig_pos = deepcopy(at.atom_data.position[i_at])
    pos = single_atom_random_walk!(pos, step_size)
    new_energy = interaction_energy(at, lj; frozen=frozen) + e_shift
    ΔE = new_energy - orig_energy
    if ΔE <= 0 && e_demon < demon_energy_threshold
        e_demon -= ΔE
        accept = true
        at = @set at.system_data.energy = new_energy
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon gains energy, accept"
    elseif 0 < ΔE <= e_demon
        e_demon -= ΔE
        accept = true
        at = @set at.system_data.energy = new_energy
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon gives energy, accept"
    else
        at = @set at.atom_data.position[i_at] = orig_pos
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon has no enough energy, reject"
    end
    return accept, at, e_demon
end

function additional_demon_walk!(
                    e_demon::Float64,
                    at::Atoms, 
                    lj::LJParameters, 
                    step_size::Float64; 
                    frozen::Int=0, 
                    e_shift::Float64=0.0, 
                    e_demon_tolerance::Float64=1e-10,
                    max_add_steps::Int=1_000_000
                    )
    accept_this_walker = false
    new_accept_count::Int = 0
    additional_steps::Int = 1
    while e_demon > e_demon_tolerance && additional_steps < max_add_steps
        accept, at, e_demon = single_atom_demon_walk!(at, lj, step_size; frozen=frozen, e_shift=e_shift, e_demon=e_demon, demon_energy_threshold=0.0)
        if accept
            new_accept_count += 1
        end
        additional_steps += 1
        if additional_steps > max_add_steps ÷ 2
            step_size /= 2
            @info "Further reduced step size to $(step_size/10)."
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
                    at::Atoms, 
                    lj::LJParameters, 
                    step_size::Float64; 
                    frozen::Int=0, 
                    e_shift::Float64=0.0, 
                    e_demon_tolerance::Float64=1e-10,
                    demon_energy_threshold::Float64=0.0,
                    max_add_steps::Int=1_000_000)
    accept_this_walker = false
    e_demon::Float64 = 0.0
    accept_count::Int64 = 0
    initial_energy::Float64  = at.system_data.energy
    @info "Initial energy: $(at.system_data.energy)"
    for i in 1:n_steps
        accept, at, e_demon = single_atom_demon_walk!(at, lj, step_size; 
                                frozen=frozen, e_shift=e_shift, e_demon=e_demon, 
                                demon_energy_threshold=demon_energy_threshold,
                                )
        if accept
            accept_count += 1
        end
    end
    if e_demon > e_demon_tolerance
        @info "Demon has too much energy remain (>tolerance=$(e_demon_tolerance)): $e_demon, performing more demon walk.
        Current accept rate: $(accept_count/n_steps). Step size is lowered."
        step_size /= 10
        accept_this_walker, _, at = additional_demon_walk!(e_demon, at, lj, step_size; 
                                                frozen=frozen, e_shift=e_shift, e_demon_tolerance=e_demon_tolerance, 
                                                max_add_steps=max_add_steps)
    else
        accept_this_walker = true
    end

    final_energy = at.system_data.energy + e_demon
    @info "Final energy: $(at.system_data.energy)"
    @debug "initial_energy: ", initial_energy, " final_energy: ", final_energy
    return accept_this_walker, accept_count/n_steps, at
end


end # module AtomsMonteCarloMoves