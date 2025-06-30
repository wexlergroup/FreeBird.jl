# Atomistic walkers
abstract type AtomWalkers <: AbstractLiveSet end

"""
    assign_energy!(walker::AtomWalker, lj::LennardJonesParametersSets)

Assigns the energy to the given `walker` using the Lennard-Jones parameters `lj`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `lj::LennardJonesParametersSets`: The Lennard-Jones parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_energy!(walker::AtomWalker, lj::LennardJonesParametersSets)
    # walker.energy_frozen_part = frozen_energy(walker.configuration, lj, walker.list_num_par, walker.frozen)
    walker.energy = interacting_energy(walker.configuration, lj, walker.list_num_par, walker.frozen) + walker.energy_frozen_part
    return walker
end

"""
    assign_frozen_energy!(walker::AtomWalker, lj::LennardJonesParametersSets)

Assigns the frozen energy to the given `walker` using the Lennard-Jones parameters `lj`.

# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `lj::LennardJonesParametersSets`: The Lennard-Jones parameters.

# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.

"""
function assign_frozen_energy!(walker::AtomWalker, lj::LennardJonesParametersSets)
    walker.energy_frozen_part = frozen_energy(walker.configuration, lj, walker.list_num_par, walker.frozen)
    return walker
end

"""
    assign_energy!(walker::AtomWalker, lj::LennardJonesParametersSets, surface::AtomWalker)
Assigns the energy to the given `walker` using the Lennard-Jones parameters `lj` with an external surface.
# Arguments
- `walker::AtomWalker`: The walker object to assign the energy to.
- `lj::LennardJonesParametersSets`: The Lennard-Jones parameters.
- `surface::AtomWalker`: The surface walker object to consider in the energy calculation.
# Returns
- `walker::AtomWalker`: The walker object with the assigned energy.
"""
function assign_energy!(walker::AtomWalker, lj::LennardJonesParametersSets, surface::AtomWalker)
    walker.energy =  interacting_energy(walker.configuration, lj, walker.list_num_par, walker.frozen, surface.configuration) + walker.energy_frozen_part
    return walker
end


"""
    struct LJAtomWalkers <: AtomWalkers

The `LJAtomWalkers` struct represents a collection of atom walkers that interact with each other using the Lennard-Jones potential.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `lj_potential::LennardJonesParametersSets`: The Lennard-Jones potential parameters. See `LennardJonesParametersSets`.

# Constructor
- `LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::LennardJonesParametersSets; assign_energy=true)`: 
    Constructs a new `LJAtomWalkers` object with the given walkers and Lennard-Jones potential parameters. If `assign_energy=true`,
    the energy of each walker is assigned using the Lennard-Jones potential.

"""
struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    lj_potential::LennardJonesParametersSets
    function LJAtomWalkers(walkers::Vector{AtomWalker{C}}, lj_potential::LennardJonesParametersSets; assign_energy=true, const_frozen_part=true) where C
        if const_frozen_part && !isempty(walkers)
            frozen_part_energy = frozen_energy(walkers[1].configuration, lj_potential, walkers[1].list_num_par, walkers[1].frozen)
        end
        if assign_energy
            Threads.@threads for walker in walkers
                if const_frozen_part
                    walker.energy_frozen_part = frozen_part_energy
                else
                    assign_frozen_energy!(walker, lj_potential)
                end
                assign_energy!(walker, lj_potential) # comes after assign_frozen_energy!
            end
        end
        return new(walkers, lj_potential)
    end
end


"""
    struct LJSurfaceWalkers <: AtomWalkers
The `LJSurfaceWalkers` struct represents a collection of atom walkers interacting through a Lennard-Jones potential, 
with the presence of an external surface object wrapped in an `AtomWalker`.

# Fields
- `walkers::Vector{AtomWalker{C}}`: A vector of atom walkers, where `C` is the number of components.
- `lj_potential::LennardJonesParametersSets`: The Lennard-Jones potential parameters.
- `surface::AtomWalker{CS}`: An atom walker representing the surface, where `CS` is the number of components of the surface.

# Constructor
- `LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                    lj_potential::LennardJonesParametersSets, 
                    surface::AtomWalker{CS}; assign_energy=true)`

    Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones potential parameters, and a single surface walker. 
    If `assign_energy=true`, the energy of each walker is assigned using the Lennard-Jones potential and the surface.

- `LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                            lj_potential::LennardJonesParametersSets, 
                            surface::AtomWalker{CS}, 
                            assign_energy_parallel::Symbol,
                            ) where C where CS`

    Constructs a new `LJSurfaceWalkers` object with the given walkers, Lennard-Jones potential parameters, and a single surface walker.
    The `assign_energy_parallel` argument determines whether to assign energy in parallel using threads (`:threads`) or distributed 
    processes (`:distributed`).
"""
struct LJSurfaceWalkers <: AtomWalkers
    walkers::Vector{AtomWalker{C}} where C
    lj_potential::LennardJonesParametersSets
    surface::AtomWalker{CS} where CS
    function LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                                lj_potential::LennardJonesParametersSets, 
                                surface::AtomWalker{CS}; 
                                assign_energy = true,
                                ) where C where CS
        update_walker!(surface, :frozen, ones(Bool, length(surface.list_num_par)))
        frozen_part_energy = surface.energy_frozen_part # assuming (only) the surface is frozen
        if assign_energy
            Threads.@threads for walker in walkers
                walker.energy_frozen_part = frozen_part_energy
                assign_energy!(walker, lj_potential, surface)
            end
        end
        return new(walkers, lj_potential, surface)
    end
end

function LJSurfaceWalkers(walkers::Vector{AtomWalker{C}}, 
                            lj_potential::LennardJonesParametersSets, 
                            surface::AtomWalker{CS}, 
                            assign_energy_parallel::Symbol,
                            ) where C where CS
    update_walker!(surface, :frozen, ones(Bool, length(surface.list_num_par)))
    frozen_part_energy = surface.energy_frozen_part
    if assign_energy_parallel == :threads
        @info "Assigning energy to walkers in parallel using $(Threads.nthreads()) threads..."
        Threads.@threads for walker in walkers
            walker.energy_frozen_part = frozen_part_energy
            assign_energy!(walker, lj_potential, surface)
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
                    assign_energy!(walker, lj_potential, surface)
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
            assign_energy!(walker, lj_potential, surface)
        end
    else
        error("Invalid parallelization option: $assign_energy_parallel. Use :threads or :distributed.")
    end
    return LJSurfaceWalkers(walkers, lj_potential, surface)
end
