"""
    FreeBirdIO

Module for input/output operations in the FreeBird package.
"""
module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV
using AtomsBase, Unitful

using ..AbstractWalkers

export read_single_config, read_configs, read_single_walker, read_walkers
export write_single_walker, write_walkers

export DataSavingStrategy, SaveEveryN
export write_df, write_df_every_n, write_walker_every_n, write_ls_every_n

export generate_initial_configs
export convert_system_to_walker, convert_walker_to_system

"""
    set_pbc(at::Atoms, pbc::Vector)

Set the periodic boundary conditions for a system of atoms.

# Arguments
- `at::Atoms`: The system of atoms.
- `pbc::Vector`: A vector of length 3 specifying the periodic boundary conditions for each dimension. Each element can be either `true` for periodic boundary conditions or `false` for Dirichlet zero boundary conditions.

# Returns
- `FlexibleSystem`: A flexible system with the specified boundary conditions.

"""
function set_pbc(at::Atoms, pbc::Vector)
    pbc_conditions = []
    for i in 1:3
        if pbc[i] == false
            push!(pbc_conditions, DirichletZero())
        elseif pbc[i] == true
            push!(pbc_conditions, Periodic())
        else
            error("Unsupported boundary condition: $(pbc[i])")
        end
    end
    pbc_conditions = Vector{BoundaryCondition}(pbc_conditions)
    return FlexibleSystem(at;boundary_conditions=pbc_conditions)
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
    num = (haskey(data, :num_frozen_part)  && resume) ? data[:num_frozen_part] : 0
    e_frozen = (haskey(data, :energy_frozen_part) && resume) ? data[:energy_frozen_part]*u"eV" : 0.0u"eV"
    at = FastSystem(at)
    return AtomWalker(at; energy=energy, iter=iter, num_frozen_part=num, energy_frozen_part=e_frozen)
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
    num = at.num_frozen_part
    e_frozen = at.energy_frozen_part.val
    return AbstractSystem(config; energy=energy, iter=iter, num_frozen_part=num, energy_frozen_part=e_frozen)
end

"""
    append_walker(filename::String, at::AtomWalker)

Append an `AtomWalker` object to a file.
"""
function append_walker(filename::String, at::AtomWalker)
    ats = read_walkers(filename)
    push!(ats, at)
    write_walkers(filename, ats)
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
function write_walkers(filename::String, ats::Vector{AtomWalker})
    flex = convert_walker_to_system.(ats)
    save_trajectory(filename::String, flex)
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
    boundary_conditions = [DirichletZero(), DirichletZero(), DirichletZero()]
    list_of_atoms = [particle_type => [rand(), rand(), rand()] for _ in 1:num_particle]
    system = periodic_system(list_of_atoms, box, fractional=true)
    flex = FlexibleSystem(system; boundary_conditions=boundary_conditions)
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
    abstract type DataSavingStrategy

Abstract type representing a strategy for saving data.
"""
abstract type DataSavingStrategy end

"""
    struct SaveEveryN <: DataSavingStrategy

SaveEveryN is a concrete subtype of DataSavingStrategy that specifies saving data every N steps.

# Fields
- `df_filename::String`: The name of the file to save the DataFrame to.
- `wk_filename::String`: The name of the file to save the atom walker to.
- `ls_filename::String`: The name of the file to save the liveset to.
- `n::Int`: The number of steps between each save.

"""
@kwdef struct SaveEveryN <: DataSavingStrategy
    df_filename::String = "output_df.csv"
    wk_filename::String = "output.traj.extxyz"
    ls_filename::String = "output.ls.extxyz"
    n::Int = 100
end

"""
    write_df(filename::String, df::DataFrame)

Write a DataFrame to a CSV file.

# Arguments
- `filename::String`: The name of the file to write to.
- `df::DataFrame`: The DataFrame to write.

"""
write_df(filename::String, df::DataFrame) = CSV.write(filename, df; append=false)

"""
    write_df_every_n(df::DataFrame, step::Int, d_strategy::SaveEveryN)

Write the DataFrame `df` to a file specified by `d_strategy.filename` every `d_strategy.n` steps.

# Arguments
- `df::DataFrame`: The DataFrame to be written.
- `step::Int`: The current step number.
- `d_strategy::SaveEveryN`: The save strategy specifying the filename and the step interval.

"""
function write_df_every_n(df::DataFrame, step::Int, d_strategy::SaveEveryN)
    if step % d_strategy.n == 0
        write_df(d_strategy.df_filename, df)
    end
end

"""
    write_walker_every_n(at::AtomWalker, step::Int, d_strategy::SaveEveryN)

Write the atom walker `at` to a file specified by `d_strategy.wk_filename` every `d_strategy.n` steps.

# Arguments
- `at::AtomWalker`: The atom walker to be written.
- `step::Int`: The current step number.
- `d_strategy::SaveEveryN`: The save strategy specifying the file name and the interval.

"""
function write_walker_every_n(at::AtomWalker, step::Int, d_strategy::SaveEveryN)
    if step % d_strategy.n == 0
        write_single_walker(d_strategy.wk_filename, at, true)
    end
end

"""
    write_ls_every_n(ls::AtomWalkers, step::Int, d_strategy::SaveEveryN)

Write the liveset `ls` to file every `n` steps, as specified by the `d_strategy`.

# Arguments
- `ls::AtomWalkers`: The liveset to be written.
- `step::Int`: The current step number.
- `d_strategy::SaveEveryN`: The save strategy specifying the frequency of writing.

"""
function write_ls_every_n(ls::AtomWalkers, step::Int, d_strategy::SaveEveryN)
    if step % d_strategy.n == 0
        write_walkers(d_strategy.ls_filename, ls.walkers)
    end
end

end # module