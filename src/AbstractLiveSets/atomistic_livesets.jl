# Atomistic walkers
abstract type AtomWalkers <: AbstractLiveSet end

"""
    assign_energy!(walker::AtomWalker, pot::AbstractPotential)

Assigns the energy to the given `walker` using the an `AbstractPotential` `pot`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `pot::AbstractPotential`: The potential parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_energy!(walker::AtomWalker, pot::AbstractPotential)
    walker.energy = interacting_energy(walker.configuration, pot, walker.list_num_par, walker.frozen) + walker.energy_frozen_part
    return walker
end

"""
    assign_frozen_energy!(walker::AtomWalker, pot::AbstractPotential)

Assigns the frozen energy to the given `walker` using the an `AbstractPotential` `pot`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `pot::AbstractPotential`: The potential parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_frozen_energy!(walker::AtomWalker, pot::AbstractPotential)
    walker.energy_frozen_part = frozen_energy(walker.configuration, pot, walker.list_num_par, walker.frozen)
    return walker
end

"""
    assign_energy!(walker::AtomWalker, pot::AbstractPotential, surface::AtomWalker)

Assigns the energy to the given `walker` using the an `AbstractPotential` `pot` with an external surface.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `pot::AbstractPotential`: The potential parameters.
- `surface::AtomWalker`: The surface walker object to consider in the energy calculation.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.
"""
function assign_energy!(walker::AtomWalker, pot::AbstractPotential, surface::AtomWalker)
    walker.energy =  interacting_energy(walker.configuration, pot, walker.list_num_par, walker.frozen, surface.configuration) + walker.energy_frozen_part
    return walker
end

"""
    assign_energy!(walker::Vector{AtomWalker{C}}, pot::AbstractPotential) where C

Assigns the energy to each walker in `walker` using an `AbstractPotential` `pot`.

# Arguments
- `walker::Vector{AtomWalker{C}}`: A vector of walker objects to assign the energy to, where `C` is the number of components.
- `pot::AbstractPotential`: The abstract potential to use for energy assignment.

# Returns
- `walker::Vector{AtomWalker{C}}`: The vector of walker objects with the assigned energy.
"""
function assign_energy!(walker::Vector{AtomWalker{C}}, pot::AbstractPotential) where C
    for w in walker
        w.energy = interacting_energy(w.configuration, pot)
    end
    return walker
end

"""
    assign_energy!(walkers::Vector{AtomWalker{C}}, pot::Union{LJParameters, CompositeParameterSets{C, LJParameters}}; assign_energy=true, const_frozen_part=true) where C

Assigns the energy to each walker in `walkers` using the a single-component or multi-component Lennard-Jones potential `pot`.
If `const_frozen_part=true`, the frozen part of the energy is calculated only once for the first walker and assigned to all walkers.
If `assign_energy=true`, the energy is assigned to each walker.

# Arguments
- `walkers::Vector{AtomWalker{C}}`: A vector of walker objects to assign the energy to, where `C` is the number of components.
- `pot::Union{LJParameters, CompositeParameterSets{C, LJParameters}}`: The potential parameters.
- `assign_energy::Bool=true`: Whether to assign the energy to each walker.
- `const_frozen_part::Bool=true`: Whether to calculate the frozen part of the energy only once for the first walker and assign it to all walkers.

# Returns
- `walkers::Vector{AtomWalker{C}}`: The vector of walker objects with the assigned energy.
"""
function assign_energy!(walkers::Vector{AtomWalker{C}}, 
                        pot::Union{LJParameters, CompositeParameterSets{C, LJParameters}}; 
                        assign_energy=true, 
                        const_frozen_part=true
                        ) where C
    if const_frozen_part && !isempty(walkers)
        frozen_part_energy = frozen_energy(walkers[1].configuration, pot, walkers[1].list_num_par, walkers[1].frozen)
    end
    if assign_energy
        Threads.@threads for walker in walkers
            if const_frozen_part
                walker.energy_frozen_part = frozen_part_energy
            else
                assign_frozen_energy!(walker, pot)
            end
            assign_energy!(walker, pot) # comes after assign_frozen_energy!
        end
    end
    return walkers
end


"""
    struct LJAtomWalkers <: AtomWalkers

The `LJAtomWalkers` struct represents a collection of atom walkers that interact with each other using the Lennard-Jones potential.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `potential::AbstractPotential`: The Lennard-Jones potential parameters. See `AbstractPotential`.

# Constructor
- `LJAtomWalkers(walkers::Vector{AtomWalker{C}}, pot::AbstractPotential; assign_energy=true)`: 
    Constructs a new `LJAtomWalkers` object with the given walkers and Lennard-Jones potential parameters. If `assign_energy=true`,
    the energy of each walker is assigned using the Lennard-Jones potential.

"""
struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    potential::Union{LJParameters, CompositeParameterSets{C, LJParameters}} where C
    function LJAtomWalkers(walkers::Vector{AtomWalker{C}}, pot::Union{LJParameters, CompositeParameterSets{C, LJParameters}}; assign_energy=true, const_frozen_part=true) where C
        assign_energy!(walkers, pot; assign_energy=assign_energy, const_frozen_part=const_frozen_part)
        return new(walkers, pot)
    end
end

AtomWalkers(walkers::Vector{AtomWalker{C}}, potential::Union{LJParameters, CompositeParameterSets{C, LJParameters}}; assign_energy=true, const_frozen_part=true) where C = LJAtomWalkers(walkers, potential; assign_energy=assign_energy, const_frozen_part=const_frozen_part)

"""
    struct MLIPAtomWalkers <: AtomWalkers
The `MLIPAtomWalkers` struct represents a collection of atom walkers that interact with each other using a machine learning interatomic potential (MLIP).

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `potential::PyMLPotential`: The machine learning interatomic potential wrapped in a `PyMLPotential`.

# Constructor
- `MLIPAtomWalkers(walkers::Vector{AtomWalker{C}}, pot::PyMLPotential; assign_energy=true)`: 
    Constructs a new `MLIPAtomWalkers` object with the given walkers and MLIP potential. If `assign_energy=true`,
    the energy of each walker is assigned using the MLIP potential.
"""
struct MLIPAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    potential::PyMLPotential
    function MLIPAtomWalkers(walkers::Vector{AtomWalker{C}}, pot::PyMLPotential; assign_energy=true) where {C}
        if assign_energy
            assign_energy!(walkers, pot)
        end
        return new(walkers, pot)
    end
end

AtomWalkers(walkers::Vector{AtomWalker{C}}, potential::PyMLPotential; assign_energy=true) where C = MLIPAtomWalkers(walkers, potential; assign_energy=assign_energy)


"""
    struct LJSurfaceWalkers <: AtomWalkers

The `LJSurfaceWalkers` struct represents a collection of atom walkers interacting through a Lennard-Jones potential, 
with the presence of an external surface object wrapped in an `AtomWalker`.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `pot::CompositeParameterSets{CP, LJParameters}`: The Lennard-Jones parameter sets. The surface
  energy paths evaluate walker components against the surface as an appended LAST component, so
  the parameter set must carry `CP = C + 1` components (walker components first, surface last);
  the constructors validate this and throw an `ArgumentError` on a mismatch. A bare
  `LJParameters` is not accepted: no surface energy path can evaluate it.
- `surface::AtomWalker{CS}`: An atom walker representing the surface, where `CS` is the number of components of the surface.

# The frozen-part energy convention
The constructors READ `surface.energy_frozen_part` and copy it into every walker; they do not
compute it unless asked. Callers must either assign the adsorbent self-energy before
construction (idiom: `surface.energy_frozen_part = interacting_energy(surface.configuration, lj)`)
or pass `compute_frozen_energy=true`, which computes it from the surface configuration with the
parameter set's surface-surface entry. Under the field's `0.0u"eV"` default, walker energies
silently omit the entire adsorbent self-energy: within one liveset that is a constant shift, but
any quantity comparing totals across systems or conventions inherits it as an error.

# Constructor
- `LJSurfaceWalkers(walkers::Vector{AtomWalker{C}},
                    pot::CompositeParameterSets{CP, LJParameters},
                    surface::AtomWalker{CS}; assign_energy=true, compute_frozen_energy=false)`

    Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones parameter
    sets, and a single surface walker. If `assign_energy=true`, the energy of each walker is
    assigned using the Lennard-Jones parameters and the surface.

- `LJSurfaceWalkers(walkers::Vector{AtomWalker{C}},
                            pot::CompositeParameterSets{CP, LJParameters},
                            surface::AtomWalker{CS},
                            assign_energy_parallel::Symbol;
                            compute_frozen_energy=false,
                            ) where {C, CP, CS}

    Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones parameter
    sets, and a single surface walker.
    The `assign_energy_parallel` argument determines whether to assign energy in parallel using threads (`:threads`) or distributed
    processes (`:distributed`).
"""
struct LJSurfaceWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    potential::CompositeParameterSets{CP, LJParameters} where CP
    surface::AtomWalker{CS} where CS
    function LJSurfaceWalkers(walkers::Vector{AtomWalker{C}},
                                pot::CompositeParameterSets{CP, LJParameters},
                                surface::AtomWalker{CS};
                                assign_energy = true,
                                compute_frozen_energy = false,
                                ) where {C, CP, CS}
        # The surface energy paths append the surface as the LAST component, so
        # the parameter set must carry exactly one more component than the
        # walkers (walker components first, surface last); a mismatch would
        # make the total-energy and single-site paths read different
        # parameter-set entries with no error.
        CP == C + 1 || throw(ArgumentError(
            "LJSurfaceWalkers requires a parameter set with one component per " *
            "walker component plus one for the surface (walker components " *
            "first, surface last): got $CP-component parameters over " *
            "$C-component walkers, expected $(C + 1)"))
        update_walker!(surface, :frozen, ones(Bool, length(surface.list_num_par)))
        if compute_frozen_energy
            surface.energy_frozen_part = interacting_energy(surface.configuration, pot.param_sets[end, end])
        end
        frozen_part_energy = surface.energy_frozen_part # assuming (only) the surface is frozen
        if assign_energy
            Threads.@threads for walker in walkers
                walker.energy_frozen_part = frozen_part_energy
                assign_energy!(walker, pot, surface)
            end
        end
        return new(walkers, pot, surface)
    end
end

"""
    LJSurfaceWalkers(walkers::Vector{AtomWalker{C}},
                            pot::CompositeParameterSets{CP, LJParameters},
                            surface::AtomWalker{CS},
                            assign_energy_parallel::Symbol;
                            compute_frozen_energy=false,
                            ) where {C, CP, CS}

Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones parameter sets, and a single surface walker.
The `assign_energy_parallel` argument determines whether to assign energy in parallel using threads (`:threads`) or distributed
processes (`:distributed`). The parameter set must carry `CP = C + 1` components (walker
components first, surface last; validated with an `ArgumentError`), and
`surface.energy_frozen_part` is read, not computed, unless `compute_frozen_energy=true`; see
[`LJSurfaceWalkers`](@ref).

# Arguments
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `pot::CompositeParameterSets{CP, LJParameters}`: The Lennard-Jones parameter sets.
- `surface::AtomWalker{CS}`: An atom walker representing the surface, where `CS` is the number of components of the surface.
- `assign_energy_parallel::Symbol`: The method to use for parallel energy assignment. Can be `:threads` or `:distributed`.
- `compute_frozen_energy::Bool`: Whether to compute `surface.energy_frozen_part` from the surface configuration before assigning walker energies. Default is `false` (the field is read as-is).

# Returns
- `LJSurfaceWalkers`: A new `LJSurfaceWalkers` object with the assigned energy.
"""
function LJSurfaceWalkers(walkers::Vector{AtomWalker{C}},
                            pot::CompositeParameterSets{CP, LJParameters},
                            surface::AtomWalker{CS},
                            assign_energy_parallel::Symbol;
                            compute_frozen_energy = false,
                            ) where {C, CP, CS}
    CP == C + 1 || throw(ArgumentError(
        "LJSurfaceWalkers requires a parameter set with one component per " *
        "walker component plus one for the surface (walker components " *
        "first, surface last): got $CP-component parameters over " *
        "$C-component walkers, expected $(C + 1)"))
    update_walker!(surface, :frozen, ones(Bool, length(surface.list_num_par)))
    if compute_frozen_energy
        surface.energy_frozen_part = interacting_energy(surface.configuration, pot.param_sets[end, end])
    end
    frozen_part_energy = surface.energy_frozen_part
    if assign_energy_parallel == :threads
        @info "Assigning energy to walkers in parallel using $(Threads.nthreads()) threads..."
        Threads.@threads for walker in walkers
            walker.energy_frozen_part = frozen_part_energy
            assign_energy!(walker, pot, surface)
        end
    elseif assign_energy_parallel == :distributed
        @info "Assigning energy to walkers in parallel using $(nworkers()) distributed processes..."
        current_first_task = 1
        remaining_tasks = length(walkers)
        while current_first_task  + nworkers() - 1 <= length(walkers) && remaining_tasks >= nworkers()
            spawned = Vector{Future}(undef, nworkers())
            for i in current_first_task:current_first_task + nworkers() - 1
                worker_id = workers()[mod1(i, length(workers()))]
                walker = walkers[i]
                spawned[mod1(i, length(workers()))] = @spawnat worker_id begin
                    walker.energy_frozen_part = frozen_part_energy
                    assign_energy!(walker, pot, surface)
                end
            end
            fetch.(spawned) # Wait for all workers to finish
            remaining_tasks = length(walkers) - (current_first_task + nworkers() - 1)
            current_first_task += nworkers()
            @info "remaining tasks: $remaining_tasks"
        end
        @info "Assigning energy to the remaining $(remaining_tasks) walkers..."
        for walker in walkers[end-remaining_tasks+1:end]
            walker.energy_frozen_part = frozen_part_energy
            assign_energy!(walker, pot, surface)
        end
    else
        error("Invalid parallelization option: $assign_energy_parallel. Use :threads or :distributed.")
    end
    return LJSurfaceWalkers(walkers, pot, surface)
end

"""
    struct GuptaAtomWalkers <: AtomWalkers

The `GuptaAtomWalkers` struct represents a collection of atom walkers that interact with each other using the Gupta potential.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `potential::Union{GuptaParameters, CompositeParameterSets{C, GuptaParameters}}`: The Gupta potential parameters. See `GuptaParameters` for more details.

"""
struct GuptaAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    potential::Union{GuptaParameters, CompositeParameterSets{C, GuptaParameters}} where C
    function GuptaAtomWalkers(walkers::Vector{AtomWalker{C}}, gupta_potential::Union{GuptaParameters, CompositeParameterSets{C, GuptaParameters}}; assign_energy=true) where C
        if assign_energy
            assign_energy!(walkers, gupta_potential)
        end
        return new(walkers, gupta_potential)
    end
end

AtomWalkers(walkers::Vector{AtomWalker{C}}, potential::Union{GuptaParameters, CompositeParameterSets{C, GuptaParameters}}; assign_energy=true) where C = GuptaAtomWalkers(walkers, potential; assign_energy=assign_energy)