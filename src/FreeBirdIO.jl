module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV

export read_walkers, write_walkers
export read_single_walker, write_single_walker
export write_df

read_walkers(filename::String) = Atoms.(read_frames(filename::String))

write_walkers(filename::String, ats::Vector) = save_trajectory(filename::String, ats) 

read_single_walker(filename::String) = Atoms(read_frame(filename::String))

write_single_walker(filename::String, at::Atoms) = save_system(filename::String, at)

write_df(filename::String, df::DataFrame) = CSV.write(filename, df)


end # module