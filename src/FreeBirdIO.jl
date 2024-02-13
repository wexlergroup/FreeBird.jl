module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV
using AtomsBase, Unitful

using ..AbstractWalkers

export read_single_config, read_configs, read_single_walker, read_walkers
export write_df
export generate_initial_configs
export convert_system_to_walker, convert_walker_to_system

# read_configs(filename::String) = FlexibleSystem.(Atoms.(read_frames(filename::String)))

# read_single_config(filename::String) = FlexibleSystem(Atoms(read_frame(filename::String)))

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
    pbc_conditions = Vector{BoundaryCondition}(pbc_conditions)
    return FlexibleSystem(at;boundary_conditions=pbc_conditions)
end

function convert_system_to_walker(at::FlexibleSystem)
    data = at.data
    energy = haskey(data, :energy) ? data[:energy]*u"eV" : 0.0u"eV"
    iter = haskey(data, :iter) ? data[:iter] : 0
    num = haskey(data, :num_frozen_part) ? data[:num_frozen_part] : 0
    e_frozen = haskey(data, :energy_frozen_part) ? data[:energy_frozen_part]*u"eV" : 0.0u"eV"
    at = FastSystem(at)
    return AtomWalker(at; energy=energy, iter=iter, num_frozen_part=num, energy_frozen_part=e_frozen)
end

function read_single_config(filename::String, pbc::Vector)
    at = Atoms(read_frame(filename::String))
    return set_pbc(at, pbc)
end

function read_single_config(filename::String; pbc::String="TTT")
    pbc = [pbc[1] == 'T', pbc[2] == 'T', pbc[3] == 'T']
    read_single_config(filename, pbc)
end

function read_configs(filename::String, pbc::Vector)
    ats = Atoms.(read_frames(filename::String))
    return [set_pbc(at, pbc) for at in ats]
end

function read_configs(filename::String; pbc::String="TTT")
    pbc = [pbc[1] == 'T', pbc[2] == 'T', pbc[3] == 'T']
    read_configs(filename, pbc)
end

function read_single_walker(filename::String; pbc::String="TTT")
    at = read_single_config(filename; pbc=pbc)
    return convert_system_to_walker(at)    
end

function read_walkers(filename::String; pbc::String="TTT")
    ats = read_configs(filename; pbc=pbc)
    return convert_system_to_walker.(ats)
end

function convert_walker_to_system(at::AtomWalker)
    config::FastSystem = at.configuration
    energy = at.energy.val
    iter = at.iter
    num = at.num_frozen_part
    e_frozen = at.energy_frozen_part.val
    return AbstractSystem(config; energy=energy, iter=iter, num_frozen_part=num, energy_frozen_part=e_frozen)
end

function write_single_walker(filename::String, at::AtomWalker)
    flex = convert_walker_to_system(at)
    save_system(filename::String, flex)
end

function write_walkers(filename::String, ats::Vector{AtomWalker})
    flex = convert_walker_to_system.(ats)
    save_trajectory(filename::String, flex)
end


function generate_random_starting_config(volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    # generate random starting configuration
    total_volume = volume_per_particle * num_particle
    box_length = total_volume^(1/3)
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    list_of_atoms = [particle_type => [rand(), rand(), rand()] for _ in 1:num_particle]
    system = periodic_system(list_of_atoms, box, fractional=true)
    return FastSystem(system)
end

function generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    [generate_random_starting_config(volume_per_particle, num_particle; particle_type=particle_type) for _ in 1:num_walkers]
end

write_df(filename::String, df::DataFrame) = CSV.write(filename, df; append=true)

end # module