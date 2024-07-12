"""
Module containing functions for performing Monte Carlo moves on atomic/molecular systems.
"""
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
export build_new_box


"""
    periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T

Wrap the position vector `pos` according to the periodic boundary conditions of the `system`. If the
boundary condition is `Periodic()`, the position is wrapped using the modulo operator. If the boundary
condition is `DirichletZero()`, the position is wrapped by reflecting the position vector across the boundary.

# Arguments
- `pos::SVector{3,T}`: The position vector to be wrapped.
- `system::AbstractSystem`: The system containing the periodic boundary conditions.

# Returns
The wrapped position vector.

"""
function periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T
    pbc = system.boundary_conditions
    box = system.bounding_box
    new_pos = Vector{typeof(0.0u"Å")}(undef, 3)
    for i in eachindex(pbc)
        if pbc[i] == Periodic() # wrap the position
            new_pos[i] = mod(pos[i], box[i][i])
        elseif pbc[i] == DirichletZero() # reflect the position
            if pos[i] > box[i][i]
                new_pos[i] = box[i][i]*2 - pos[i]
            elseif pos[i] < 0.0u"Å"
                new_pos[i] = -pos[i]
            else
                new_pos[i] = pos[i]
            end
        else
            error("Unsupported boundary condition: $(pbc[i]).")
        end
    end
    pos = SVector{3,T}(new_pos)
    return pos
end

function build_new_box(old_box, slab_height)
    new_box = [[old_box[i][j] for j in 1:3] for i in 1:3]
    new_box[3][3] -= slab_height
    return new_box
end


function periodic_boundary_wrap!(pos::SVector{3,T}, pbc, box, slab_h) where T
    new_pos = Vector{typeof(0.0u"Å")}(undef, 3)
    for i in eachindex(pbc)
        if pbc[i] == Periodic() # wrap the position
            if i ==3
                pos[i] = pbc[i] - slab_h
                new_pos[i] = mod(pos[i], box[i][i]) + slab_h
            elseif i != 3
                pos[i] = pbc[i]
                new_pos[i] = mod(pos[i], box[i][i])
            end

            
        elseif pbc[i] == DirichletZero() # reflect the position
            if i ==3
                if (pos[i] - slab_h) > box[i][i]
                    new_pos[i] = box[i][i]*2 - (pos[i] - slab_h) - slab_h
                elseif (pos[i] - slab_h) < 0.0u"Å"
                    new_pos[i] = -(pos[i] - slab_h) + slab_h
                else
                    new_pos[i] = (pos[i] - slab_h) + slab_h
                end
                
            elseif i != 3
                if pos[i] > box[i][i]
                    new_pos[i] = box[i][i]*2 - pos[i]
                elseif pos[i] < 0.0u"Å"
                    new_pos[i] = -pos[i]
                else
                    new_pos[i] = pos[i]
                end
                
            end

        else
            error("Unsupported boundary condition: $(pbc[i]).")
        end
    end
    pos = SVector{3,T}(new_pos)
    return pos
end

"""
    mean_sq_displacement(at::AtomWalker, at_orig::AtomWalker)

Calculate the mean squared displacement before and after random walk(s). Note that due to the current implementation of the periodic boundary wrap, this function is not appropriate to use for calculating mean displacements in a propagation.

# Arguments
- `at::AtomWalker{C}`: The current `AtomWalker` after the random walk.
- `at_orig::AtomWalker{C}`: The original `AtomWalker` before the random walk.

# Returns
- `distsq::typeof(0.0u"Å"^2)`: The mean squared displacement of all free particles.

"""
function mean_sq_displacement(at::AtomWalker{C}, at_orig::AtomWalker{C}) where C
    distsq::typeof(0.0u"Å"^2) = 0.0u"Å"^2
    free_index = free_par_index(at)
    for i in free_index
        dist::typeof(0.0u"Å") = pbc_dist(at.configuration.position[i], at_orig.configuration.position[i], at.configuration)
        distsq += dist^2
    end
    return distsq/length(free_index)
end


"""
    single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64) where T

Perform a single atom random walk by updating the position `pos` in each direction by a random amount.
The `step_size` determines the maximum distance the atom can move in any direction.

# Arguments
- `pos::SVector{3,T}`: The current position of the atom as a 3D vector.
- `step_size::Float64`: The maximum distance the atom can move in any direction.

# Returns
- `pos`: The updated position of the atom.

"""
function single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64) where T
    (dx, dy, dz) = (rand(Uniform(-step_size,step_size)) for _ in 1:3) .* unit(T)
    pos = pos .+ (dx, dy, dz)
    return pos 
end

function free_par_index(at::AtomWalker{C}) where C
    ind_free_par = Array{Int}(undef, 0)
    comp_cut = vcat([0],cumsum(at.list_num_par))
    # @show comp_cut
    for i in 1:C
        if !at.frozen[i]
            ind_free_par = vcat(ind_free_par,comp_cut[i]+1:comp_cut[i+1])
        end
    end
    return ind_free_par
end

"""
    MC_random_walk!(n_steps::Int, at::AtomWalker, lj::LJParameters, step_size::Float64, emax::typeof(0.0u"eV"))

Perform a Monte Carlo random walk on the atomic/molecular system.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `lj::LennardJonesParametersSets`: The Lennard-Jones potential parameters.
- `step_size::Float64`: The maximum distance an atom can move in any direction.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `at::AtomWalker`: The updated walker.

"""
function MC_random_walk!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    lj::LennardJonesParametersSets, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV");
                    slab_height = 0.0u"Å",
                    dispersion_params = []
                    ) where C
    n_accept = 0
    accept_this_walker = false
    pbc = at.configuration.boundary_conditions
    box = build_new_box(at.configuration.bounding_box, slab_height)
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, pbc, box, slab_height)
        config.position[i_at] = pos
        energy = interacting_energy(config, lj, at.list_num_par, at.frozen, dispersion_params=dispersion_params) + at.energy_frozen_part
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


"""
    single_atom_demon_walk!(at::AtomWalker{C}, 
                            lj::LennardJonesParametersSets, 
                            step_size::Float64;
                            e_demon=0.0u"eV",
                            demon_energy_threshold=Inf*u"eV",
                            demon_gain_threshold=Inf*u"eV")

Perform a single atom demon walk.

Arguments:
- `at::AtomWalker`: The walker to perform the demon walk on.
- `lj::LJParameters`: The LJ parameters.
- `step_size::Float64`: The step size for the random walk.

Keyword Arguments:
- `e_demon=0.0u"eV"`: The energy of the demon.
- `demon_energy_threshold=Inf*u"eV"`: The energy threshold for the demon.
- `demon_gain_threshold=Inf*u"eV"`: The energy gain threshold for the demon during each move.

Returns:
- `accept::Bool`: Whether the move is accepted or rejected.
- `at::AtomWalker`: The updated atom walker object.
- `e_demon::Float64`: The updated energy of the demon.
"""
function single_atom_demon_walk!(
                        at::AtomWalker{C}, 
                        lj::LennardJonesParametersSets, 
                        step_size::Float64;
                        e_demon=0.0u"eV",
                        demon_energy_threshold=Inf*u"eV",
                        demon_gain_threshold=Inf*u"eV",
                        ) where C
    accept = false
    orig_energy = at.energy
    config = at.configuration
    free_index = free_par_index(at)
    i_at = rand(free_index)
    pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
    orig_pos = deepcopy(pos)
    pos = single_atom_random_walk!(pos, step_size)
    pos = periodic_boundary_wrap!(pos, config)
    config.position[i_at] = pos
    new_energy = interacting_energy(config, lj, at.list_num_par, at.frozen) + at.energy_frozen_part
    ΔE = new_energy - orig_energy
    if -demon_gain_threshold <= ΔE <= 0.0u"eV" && e_demon-ΔE < demon_energy_threshold
        e_demon += -ΔE
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
        @debug "ΔE: ", ΔE, " e_demon: ", e_demon, "demon has insufficient energy, reject"
    end
    return accept, at, e_demon
end

"""
    additional_demon_walk!(e_demon::typeof(0.0u"eV"), at::AtomWalker, lj::LJParameters, step_size::Float64;
                          e_demon_tolerance=1e-9u"eV", max_add_steps::Int=1_000_000)

Performs additional demon walk steps until the demon energy `e_demon` is below the tolerance `e_demon_tolerance`
or the maximum number of additional steps `max_add_steps` is reached.

# Arguments
- `e_demon::typeof(0.0u"eV")`: The initial demon energy.
- `at::AtomWalker{C}`: The walker that the demon walk is performed on.
- `lj::LennardJonesParametersSets`: The LJ parameters.
- `step_size::Float64`: The step size for the demon walk.
- `e_demon_tolerance=1e-9u"eV"`: The tolerance for the demon energy.
- `max_add_steps::Int=1_000_000`: The maximum number of additional steps.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted after the additional demon walk.
- `accept_rate::Float64`: The acceptance rate of the additional demon walk.
- `at::AtomWalker`: The updated walker.
- `e_demon::typeof(0.0u"eV")`: The final demon energy.

"""
function additional_demon_walk!(
                    e_demon::typeof(0.0u"eV"),
                    at::AtomWalker{C},
                    lj::LennardJonesParametersSets, 
                    step_size::Float64;
                    e_demon_tolerance=1e-9u"eV",
                    max_add_steps::Int=1_000_000
                    ) where C
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

"""
    MC_nve_walk!(
        n_steps::Int, 
        at::AtomWalker{C}, 
        lj::LennardJonesParametersSets, 
        step_size::Float64; 
        e_demon_tolerance=1e-9u"eV",
        demon_energy_threshold=Inf*u"eV",
        demon_gain_threshold=Inf*u"eV",
        max_add_steps::Int=1_000_000
    )

Perform a NVE walk using the demon algorithm. The demon algorithm is used to maintain the energy (E) of the system.


# Arguments
- `n_steps::Int`: The number of demon walks to perform.
- `at::AtomWalker`: The walker to perform the demon walk on.
- `lj::LJParameters`: The Lennard-Jones parameters for the system.
- `step_size::Float64`: The step size for the demon walk.

# Optional Arguments
- `e_demon_tolerance=1e-9u"eV"`: The energy tolerance for the demon, below which the demon walk is considered successful, 
    i.e., the NVE condition is satisfied.
- `demon_energy_threshold=Inf*u"eV"`: The maximum energy allowed for the demon to have.
- `demon_gain_threshold=Inf*u"eV"`: The energy gain threshold for the demon during each move.
- `max_add_steps::Int=1_000_000`: The maximum number of additional demon walks if the demon energy is above the tolerance.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_ratio::Float64`: The acceptance ratio of the demon walk (excluding additional demon walks).
- `at_final::AtomWalker`: The final walker after the demon walks.
- `demon_energies::Array{Float64}`: The energies of the demon at each step.
- `temp_estimate::Float64`: The estimated temperature of the system.

"""
function MC_nve_walk!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    lj::LennardJonesParametersSets, 
                    step_size::Float64; 
                    e_demon_tolerance=1e-9u"eV",
                    demon_energy_threshold=Inf*u"eV",
                    demon_gain_threshold=Inf*u"eV",
                    max_add_steps::Int=1_000_000) where C
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
    temp_per_K = temp_estimate/(8.617333262e-5u"eV/K")
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