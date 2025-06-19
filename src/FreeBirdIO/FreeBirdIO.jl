"""
    FreeBirdIO

Module for input/output operations in the FreeBird package.
"""
module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV
using AtomsBase, Unitful
using Arrow

using ..AbstractWalkers
using ..AbstractLiveSets
using ..EnergyEval

export read_single_config, read_configs, read_single_walker, read_walkers
export write_single_walker, write_walkers

export append_system

export DataSavingStrategy, SaveEveryN, SaveFreePartEveryN
export write_df, write_df_every_n, write_walker_every_n, write_ls_every_n

export generate_initial_configs
export convert_system_to_walker, convert_walker_to_system

# include the save strategies
include("save_strategies.jl")

"""
    set_pbc(at::Atoms, pbc::Vector)

Set the periodic boundary conditions for a system of atoms.

# Arguments
- `at::Atoms`: The system of atoms.
- `pbc::Vector`: A vector of length 3 specifying the periodic boundary conditions for each dimension. Each element can be either `true` for periodic boundary conditions or `false` for Dirichlet zero boundary conditions.

# Returns
- `FlexibleSystem`: A flexible system with the specified boundary conditions.

"""
function set_pbc(at::Atoms, pbc::Vector{Bool})
    return FlexibleSystem(at, periodicity=Tuple(pbc))
end

"""
    convert_system_to_walker(at::FlexibleSystem, resume::Bool)

Converts a `FlexibleSystem` object to an `AtomWalker` object.

# Arguments
- `at::FlexibleSystem`: The `FlexibleSystem` object to convert.
- `resume::Bool`: Whether to resume from previous data.

# Returns
- `AtomWalker`: The converted `AtomWalker` object.

"""
function convert_system_to_walker(at::FlexibleSystem, resume::Bool)
    data = at.data
    energy = (haskey(data, :energy) && resume) ? data[:energy]*u"eV" : 0.0u"eV"
    iter = (haskey(data, :iter) && resume && data[:iter]>=0) ? data[:iter] : 0
    list_num_par = (haskey(data, :list_num_par) && resume) ? data[:list_num_par] : [length(at)]
    C = length(list_num_par)
    frozen = (haskey(data, :frozen) && resume) ? data[:frozen] : zeros(Bool, C)
    e_frozen = (haskey(data, :energy_frozen_part) && resume) ? data[:energy_frozen_part]*u"eV" : 0.0u"eV"
    new_list = [Atom(atomic_symbol(i),position(i)) for i in at.particles]
    at = FastSystem(new_list, cell_vectors(at), periodicity(at))
    return AtomWalker{C}(at; energy=energy, iter=iter, list_num_par=list_num_par, frozen=frozen, energy_frozen_part=e_frozen)
end

"""
    read_single_config(filename::String, pbc::Vector)

Reads a single configuration from the specified file and sets the periodic boundary conditions (PBC) for the atoms.

# Arguments
- `filename::String`: The name of the file to read the configuration from.
- `pbc::Vector`: A vector specifying the periodic boundary conditions.

# Returns
- `at::Atoms`: The atoms with the PBC set.

"""
function read_single_config(filename::String, pbc::Vector)
    at = Atoms(read_frame(filename::String))
    return set_pbc(at, pbc)
end

"""
    read_single_config(filename::String; pbc::String="TTT")

Reads a single configuration from the specified file.

# Arguments
- `filename::String`: The name of the file to read from.
- `pbc::String="TTT"`: The periodic boundary conditions. Default is "TTT".

# Returns
- The configuration read from the file.

"""
function read_single_config(filename::String; pbc::String="TTT")
    pbc = [pbc[1] == 'T', pbc[2] == 'T', pbc[3] == 'T']
    read_single_config(filename, pbc)
end

"""
    read_configs(filename::String, pbc::Vector)

Reads atomic configurations from a file and applies periodic boundary conditions.

# Arguments
- `filename::String`: The name of the file containing the atomic configurations.
- `pbc::Vector`: A vector specifying the periodic boundary conditions.

# Returns
An array of atomic configurations with periodic boundary conditions applied.

"""
function read_configs(filename::String, pbc::Vector)
    ats = Atoms.(read_frames(filename::String))
    return [set_pbc(at, pbc) for at in ats]
end

"""
    read_configs(filename::String; pbc::String="TTT")

Reads configurations from a file.

# Arguments
- `filename::String`: The name of the file to read configurations from.
- `pbc::String="TTT"`: Periodic boundary conditions. A string of length 3, where each character represents whether the corresponding dimension has periodic boundary conditions ('T') or not ('F').

# Returns
- The configurations read from the file.

"""
function read_configs(filename::String; pbc::String="TTT")
    pbc = [pbc[1] == 'T', pbc[2] == 'T', pbc[3] == 'T']
    read_configs(filename, pbc)
end

"""
    read_single_walker(filename::String; pbc::String="TTT", resume::Bool=true)

Reads a single walker from the specified file.

# Arguments
- `filename::String`: The path to the file containing the walker data.
- `pbc::String`: (optional) The periodic boundary conditions. Default is "TTT".
- `resume::Bool`: (optional) Whether to resume reading from a previous checkpoint. Default is `true`.

# Returns
- The walker object read from the file.

"""
function read_single_walker(filename::String; pbc::String="TTT", resume::Bool=true)
    at = read_single_config(filename; pbc=pbc)
    return convert_system_to_walker(at, resume)    
end

"""
    read_walkers(filename::String; pbc::String="TTT", resume::Bool=true)

Reads walker configurations from a file.

# Arguments
- `filename::String`: The name of the file to read the walker configurations from.
- `pbc::String`: A string specifying the periodic boundary conditions. Default is "TTT".
- `resume::Bool`: A boolean indicating whether to resume reading from a previous checkpoint. Default is `true`.

# Returns
An array of walker objects.

"""
function read_walkers(filename::String; pbc::String="TTT", resume::Bool=true)
    ats = read_configs(filename; pbc=pbc)
    return convert_system_to_walker.(ats, resume)
end

"""
    convert_walker_to_system(at::AtomWalker)

Converts an `AtomWalker` object to an `AbstractSystem` object.

# Arguments
- `at::AtomWalker`: The `AtomWalker` object to be converted.

# Returns
- `AbstractSystem`: The converted `AbstractSystem` object.

"""
function convert_walker_to_system(at::AtomWalker)
    config::FastSystem = at.configuration
    energy = at.energy.val
    iter = at.iter
    num = join(at.list_num_par, ",")
    frozen = join(at.frozen, ",")
    e_frozen = at.energy_frozen_part.val
    return AbstractSystem(config; energy=energy, iter=iter, list_num_par=num, frozen=frozen, energy_frozen_part=e_frozen)
end

"""
    append_walker(filename::String, at::AtomWalker)

Append an `AtomWalker` object to a file.
"""
function append_walker(filename::String, at::AtomWalker{C}) where C
    sys = convert_walker_to_system(at)
    write_frame(filename, ExtXYZ.write_dict(Atoms(sys)), append=true)
end

"""
    append_system(ats1::FlexibleSystem, ats2::FlexibleSystem)

Append two `FlexibleSystem` objects into a single `FastSystem` object.
The first argument is the system to be appended to, and its bounding box and boundary conditions 
will be used for the new system.

# Arguments
- `ats1::FlexibleSystem`: The base system to be appended.
- `ats2::FlexibleSystem`: The system to append.

# Returns
- `new_list`: A new `FastSystem` object containing the appended systems.

"""
function append_system(ats1::FlexibleSystem, ats2::FlexibleSystem)
    new_list = [Atom(atomic_symbol(i),position(i)) for i in ats1.particles]
    append!(new_list,[Atom(atomic_symbol(i),position(i)) for i in ats2.particles])
    return FastSystem(new_list, cell_vectors(ats1), periodicity(ats1))
end

"""
    append_system(ats1::FlexibleSystem, ats2::Vector{FlexibleSystem})
Append a `FlexibleSystem` object to a vector of `FlexibleSystem` objects.
The first argument is the system to be appended to, and its bounding box and boundary conditions
will be used for the new systems.
# Arguments
- `ats1::FlexibleSystem`: The base system to be appended.
- `ats2::Vector{FlexibleSystem}`: A vector of `FlexibleSystem` objects to append.
# Returns
- `configs`: A vector of `FastSystem` objects containing the appended systems.
"""
function append_system(ats1::FlexibleSystem, ats2::Vector{T}) where T
    configs = Vector{FastSystem}(undef, length(ats2))
    Threads.@threads for i in eachindex(ats2)
        configs[i] = append_system(ats1, ats2[i])
    end
    return configs
end

"""
    write_single_walker(filename::String, at::AtomWalker)

Write a single AtomWalker object to a file.

# Arguments
- `filename::String`: The name of the file to write to.
- `at::AtomWalker`: The AtomWalker object to write.

"""
function write_single_walker(filename::String, at::AtomWalker)
    flex = convert_walker_to_system(at)
    save_system(filename::String, flex)
end

"""
    write_walkers(filename::String, ats::Vector{AtomWalker})

Write a collection of `AtomWalker` objects to a file.

# Arguments
- `filename::String`: The name of the file to write the walkers to.
- `ats::Vector{AtomWalker}`: The collection of `AtomWalker` objects to write.

"""
function write_walkers(filename::String, ats::Vector{AtomWalker{C}}) where C
    flex = convert_walker_to_system.(ats)
    save_trajectory(filename::String, flex)
end
"""
    write_walkers(filename::String, ats::Vector{LatticeWalker})

Write a collection of `LatticeWalker` objects to a file.

# Arguments
- `filename::String`: The name of the file to write the walkers to.
- `ats::Vector{LatticeWalker}`: The collection of `LatticeWalker` objects to write.

"""
function write_walkers(filename::String, ats::Vector{LatticeWalker{C}}) where C
    occupancies = [at.configuration.components for at in ats]
    energies = [at.energy.val for at in ats]
    energy_and_unit = "energy_$(unit(ats[1].energy))"
    df = DataFrame(energy_and_unit=>energies, "components"=>occupancies)
    CSV.write(filename, df; append=false)
end

"""
    write_single_walker(filename::String, at::AtomWalker, append::Bool)

Write a single AtomWalker object to a file. If the file already exists, append the walker to the file.

# Arguments
- `filename::String`: The name of the file to write to.
- `at::AtomWalker`: The AtomWalker object to write.
- `append::Bool`: A boolean indicating whether to append the walker to the file if it already exists.

"""
function write_single_walker(filename::String, at::AtomWalker, append::Bool)
    if isfile(filename) && append
        append_walker(filename, at)
    else
        write_walkers(filename, [at])
    end
end

"""
    extract_free_par(walker::AtomWalker)

Extract free particles from existing walker and create new walker containing only the free particles.

# Arguments
- `walker::AtomWalker`: The AtomWalker object for extraction.

"""

function extract_free_par(walker::AtomWalker)
    free_part = []
    free_indices = []
    components = split_components(walker.configuration, walker.list_num_par)
    for ind in eachindex(components)
        if !walker.frozen[ind]
            push!(free_indices, length(components[ind]))
            for i in 1:length(components[ind])
                pos = components[ind].position[i]
                sym = atomic_symbol(components[ind], i)
                push!(free_part, Atom(sym, pos))
            end
        end
    end
    fast = FastSystem(free_part, cell_vectors(walker.configuration), periodicity(walker.configuration))
    return AtomWalker{length(free_indices)}(fast;list_num_par = free_indices, energy = walker.energy - walker.energy_frozen_part, iter = walker.iter)
end

"""
    generate_random_starting_config(volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)

Generate a random starting configuration for a system of particles.

# Arguments
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Int`: The number of particles.
- `particle_type::Symbol=:H`: The type of particle (default is hydrogen).

# Returns
- `FastSystem`: A FastSystem object representing the generated system.

"""
function generate_random_starting_config(volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    # generate random starting configuration
    total_volume = volume_per_particle * num_particle
    box_length = total_volume^(1/3)
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    boundary_conditions = (false, false, false)
    list_of_atoms = [particle_type => [rand(), rand(), rand()] .* box_length * u"Å" for _ in 1:num_particle]
    system = atomic_system(list_of_atoms, box, boundary_conditions)
    return FastSystem(system)
end

"""
    generate_multi_type_random_starting_config(volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])

Generate a random starting configuration for a system of particles with multiple types.

# Arguments
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Vector{Int}`: The number of particles of each type.
- `particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)]`: The types of particles.

# Returns
- `FastSystem`: A FastSystem object representing the generated system.

"""
function generate_multi_type_random_starting_config(volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])
    num_types = length(num_particle)
    total_num_particle = sum(num_particle)
    total_volume = volume_per_particle * total_num_particle
    box_length = total_volume^(1/3)
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    boundary_conditions = (false, false, false)
    list_of_atoms = []
    for i in 1:num_types
        for _ in 1:num_particle[i]
            push!(list_of_atoms, particle_types[i] => [rand(), rand(), rand()] .* box_length * u"Å")
        end
    end
    flex = atomic_system(list_of_atoms, box, boundary_conditions)
    return FastSystem(flex)
end

"""
    generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)

Generate initial configurations for a given number of walkers.

# Arguments
- `num_walkers::Int`: The number of walkers.
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Int`: The number of particles.
- `particle_type::Symbol=:H`: The type of particle (default is :H).

# Returns
An array of initial configurations for each walker.

"""
function generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    [generate_random_starting_config(volume_per_particle, num_particle; particle_type=particle_type) for _ in 1:num_walkers]
end


"""
    generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])
Generate initial configurations for a given number of walkers with multiple particle types.
# Arguments
- `num_walkers::Int`: The number of walkers.
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Vector{Int}`: A vector containing the number of particles of each type.
- `particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)]`: A vector of symbols representing the types of particles (default is hydrogen and oxygen).
# Returns
An array of initial configurations for each walker, where each configuration contains particles of the specified types.
"""
function generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])
    [generate_multi_type_random_starting_config(volume_per_particle, num_particle; particle_types=particle_types) for _ in 1:num_walkers]
end

end # module