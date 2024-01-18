module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV

export read_liveset, save_liveset, read_walker, save_walker
export save_df

read_liveset(filename::String) = load_trajectory(filename::String)

save_liveset(filename::String, ats::Vector) = save_trajectory(filename::String, ats) 

read_walker(filename::String) = load_system(filename::String)

save_walker(filename::String, at::Atoms) = save_system(filename::String, at)

save_df(filename::String, df::DataFrame) = CSV.write(filename, df)


end # module