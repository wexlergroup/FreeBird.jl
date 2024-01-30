module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV
using AtomsBase, Unitful

export read_walkers, write_walkers
export read_single_walker, write_single_walker
export write_df
export generate_initial_liveset

read_walkers(filename::String) = FlexibleSystem.(Atoms.(read_frames(filename::String)))

write_walkers(filename::String, ats::Vector) = save_trajectory(filename::String, ats) 

read_single_walker(filename::String) = FlexibleSystem(Atoms(read_frame(filename::String)))

write_single_walker(filename::String, at::AbstractSystem) = save_system(filename::String, at)

write_df(filename::String, df::DataFrame) = CSV.write(filename, df; append=true)


function generate_random_starting_config(volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    # generate random starting configuration
    total_volume = volume_per_particle * num_particle
    box_length = total_volume^(1/3)
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    list_of_atoms = [particle_type => [rand(), rand(), rand()] for _ in 1:num_particle]
    system = periodic_system(list_of_atoms, box, fractional=true)
    return FlexibleSystem(Atoms(system))
end

function generate_initial_liveset(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    [generate_random_starting_config(volume_per_particle, num_particle; particle_type=particle_type) for _ in 1:num_walkers]
end

end # module