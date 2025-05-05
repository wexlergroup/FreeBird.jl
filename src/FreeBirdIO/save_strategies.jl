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
- `n_traj::Int`: The number of steps between each save of the culled walker into a trajectory file.
- `n_snap::Int`: The number of steps between each save of the liveset into a snapshot file.
- `n_info::Int`: The number of steps between each print of information.
"""
@kwdef struct SaveEveryN <: DataSavingStrategy
    df_filename::String = "output_df.csv"
    wk_filename::String = "output.traj.extxyz"
    ls_filename::String = "output.ls.extxyz"
    n_traj::Int = 100
    n_snap::Int = 1000
    n_info::Int = 1
end

"""
    struct SaveFreePartEveryN <: DataSavingStrategy

SaveFreePartEveryN is a concrete subtype of DataSavingStrategy that specifies saving data every N steps.
Only the free particles are saved into the trajectory and snapshot files.

# Fields
- `df_filename::String`: The name of the file to save the DataFrame to.
- `wk_filename::String`: The name of the file to save the atom walker to.
- `ls_filename::String`: The name of the file to save the liveset to.
- `n_traj::Int`: The number of steps between each save of the culled walker into a trajectory file.
- `n_snap::Int`: The number of steps between each save of the liveset into a snapshot file.
- `n_info::Int`: The number of steps between each print of information.
"""
@kwdef struct SaveFreePartEveryN <: DataSavingStrategy
    df_filename::String = "output_df.csv"
    wk_filename::String = "output.traj.extxyz"
    ls_filename::String = "output.ls.extxyz"
    n_traj::Int = 100
    n_snap::Int = 1000
    n_info::Int = 1
end

"""
    write_df(filename::String, df::DataFrame)

Write a DataFrame to a CSV/Arrow file.

# Arguments
- `filename::String`: The name of the file to write to.
- `df::DataFrame`: The DataFrame to write.

"""
function write_df(filename::String, df::DataFrame)
    if splitext(filename)[end] == ".csv"
        CSV.write(filename, df)
    elseif splitext(filename)[end] == ".arrow"
        Arrow.write(filename, df)
    else
        error("Unsupported file format. Only CSV is supported.")
    end
end

"""
    write_df_every_n(df::DataFrame, step::Int, d_strategy::SaveEveryN)

Write the DataFrame `df` to a file specified by `d_strategy.filename` every `d_strategy.n` steps.

# Arguments
- `df::DataFrame`: The DataFrame to be written.
- `step::Int`: The current step number.
- `d_strategy::SaveEveryN`: The save strategy specifying the filename and the step interval.

"""
function write_df_every_n(df::DataFrame, step::Int, d_strategy::DataSavingStrategy)
    if step % d_strategy.n_traj == 0
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
    if step % d_strategy.n_traj == 0
        write_single_walker(d_strategy.wk_filename, at, true)
    end
end

function write_walker_every_n(at::AtomWalker, step::Int, d_strategy::SaveFreePartEveryN)
    if step % d_strategy.n_traj == 0
        at = extract_free_par(at)
        write_single_walker(d_strategy.wk_filename, at, true)
    end
end

"""
    write_ls_every_n(ls::AtomWalkers, step::Int, d_strategy::SaveEveryN)

Write the liveset `ls` to file every `n` steps, as specified by the `d_strategy`.

# Arguments
- `ls::AbstractLiveSet`: The liveset to be written.
- `step::Int`: The current step number.
- `d_strategy::SaveEveryN`: The save strategy specifying the frequency of writing.

"""
function write_ls_every_n(ls::AbstractLiveSet, step::Int, d_strategy::SaveEveryN)
    if step % d_strategy.n_snap == 0
        write_walkers(d_strategy.ls_filename, ls.walkers)
    end
end

function write_ls_every_n(ls::AbstractLiveSet, step::Int, d_strategy::SaveFreePartEveryN)
    if step % d_strategy.n_snap == 0
        wks = extract_free_par.(ls.walkers)
        write_walkers(d_strategy.ls_filename, wks)
    end
end