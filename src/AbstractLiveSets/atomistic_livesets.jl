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


"""
    struct LJSurfaceWalkers <: AtomWalkers

The `LJSurfaceWalkers` struct represents a collection of atom walkers interacting through a Lennard-Jones potential, 
with the presence of an external surface object wrapped in an `AtomWalker`.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}}`: The Lennard-Jones potential parameters.
- `surface::AtomWalker{CS}`: An atom walker representing the surface, where `CS` is the number of components of the surface.

# Constructor
- `LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                    pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}}, 
                    surface::AtomWalker{CS}; assign_energy=true)`

    Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones potential parameters, and a single surface walker. 
    If `assign_energy=true`, the energy of each walker is assigned using the Lennard-Jones potential and the surface.

- `LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                            pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}},
                            surface::AtomWalker{CS}, 
                            assign_energy_parallel::Symbol,
                            ) where {C, CP, CS}

    Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones potential parameters, and a single surface walker.
    The `assign_energy_parallel` argument determines whether to assign energy in parallel using threads (`:threads`) or distributed 
    processes (`:distributed`).
"""
struct LJSurfaceWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    potential::Union{LJParameters, CompositeParameterSets{CP, LJParameters}} where CP
    surface::AtomWalker{CS} where CS
    function LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                                pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}}, 
                                surface::AtomWalker{CS}; 
                                assign_energy = true,
                                ) where {C, CP, CS}
        update_walker!(surface, :frozen, ones(Bool, length(surface.list_num_par)))
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
                            pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}},
                            surface::AtomWalker{CS}, 
                            assign_energy_parallel::Symbol,
                            ) where C where CS

Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones potential parameters, and a single surface walker.
The `assign_energy_parallel` argument determines whether to assign energy in parallel using threads (`:threads`) or distributed 
processes (`:distributed`).

# Arguments
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}}`: The Lennard-Jones potential parameters.
- `surface::AtomWalker{CS}`: An atom walker representing the surface, where `CS` is the number of components of the surface.
- `assign_energy_parallel::Symbol`: The method to use for parallel energy assignment. Can be `:threads` or `:distributed`.

# Returns
- `LJSurfaceWalkers`: A new `LJSurfaceWalkers` object with the assigned energy.
"""
function LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                            pot::Union{LJParameters, CompositeParameterSets{CP, LJParameters}}, 
                            surface::AtomWalker{CS}, 
                            assign_energy_parallel::Symbol,
                            ) where {C, CP, CS}
    update_walker!(surface, :frozen, ones(Bool, length(surface.list_num_par)))
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
