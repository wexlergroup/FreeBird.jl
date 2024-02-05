module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV
using AtomsBase, Unitful

using ..AbstractWalkers

export read_walkers, write_walkers
export read_single_walker, write_single_walker
export write_df
export generate_initial_liveset

read_walkers(filename::String) = FastSystem.(Atoms.(read_frames(filename::String)))

write_walkers(filename::String, ats::Vector) = save_trajectory(filename::String, ats) 

read_single_walker(filename::String) = FastSystem(Atoms(read_frame(filename::String)))

# write_single_walker(filename::String, at::AbstractSystem) = save_system(filename::String, at)

write_df(filename::String, df::DataFrame) = CSV.write(filename, df; append=true)

function convert_walker_to_system(at::AtomWalker)
    config::FastSystem = at.configuration
    energy = at.energy.val
    iter = at.iter
    num = at.num_frozen_part
    e_frozen = at.energy_frozen_part.val
    flex = AbstractSystem(config; energy=energy, iter=iter, num_frozen_part=num, energy_frozen_part=e_frozen)
    return flex
end

function write_single_walker(filename::String, at::AtomWalker)
    flex = convert_walker_to_system(at)
    save_system(filename::String, flex)
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

function generate_initial_liveset(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    [generate_random_starting_config(volume_per_particle, num_particle; particle_type=particle_type) for _ in 1:num_walkers]
end

end # module