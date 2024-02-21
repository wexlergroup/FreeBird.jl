var documenterSearchIndex = {"docs":
[{"location":"AtomsMCMoves/#AtomsMCMoves","page":"AtomsMCMoves","title":"AtomsMCMoves","text":"","category":"section"},{"location":"AtomsMCMoves/#Functions","page":"AtomsMCMoves","title":"Functions","text":"","category":"section"},{"location":"AtomsMCMoves/","page":"AtomsMCMoves","title":"AtomsMCMoves","text":"Modules = [AtomsMCMoves]","category":"page"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves","text":"Module containing functions for performing Monte Carlo moves on atomic/molecular systems.\n\n\n\n\n\n","category":"module"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.MC_nve_walk!-Tuple{Int64, AtomWalker, LJParameters, Float64}","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.MC_nve_walk!","text":"MC_nve_walk!(\n    n_steps::Int, \n    at::AtomWalker, \n    lj::LJParameters, \n    step_size::Float64; \n    e_demon_tolerance=1e-9u\"eV\",\n    demon_energy_threshold=Inf*u\"eV\",\n    demon_gain_threshold=Inf*u\"eV\",\n    max_add_steps::Int=1_000_000\n)\n\nPerform a NVE walk using the demon algorithm. The demon algorithm is used to maintain the energy (E) of the system.\n\nArguments\n\nn_steps::Int: The number of demon walks to perform.\nat::AtomWalker: The walker to perform the demon walk on.\nlj::LJParameters: The Lennard-Jones parameters for the system.\nstep_size::Float64: The step size for the demon walk.\n\nOptional Arguments\n\ne_demon_tolerance=1e-9u\"eV\": The energy tolerance for the demon, below which the demon walk is considered successful,    i.e., the NVE condition is satisfied.\ndemon_energy_threshold=Inf*u\"eV\": The maximum energy allowed for the demon to have.\ndemon_gain_threshold=Inf*u\"eV\": The energy gain threshold for the demon during each move.\nmax_add_steps::Int=1_000_000: The maximum number of additional demon walks if the demon energy is above the tolerance.\n\nReturns\n\naccept_this_walker::Bool: Whether the walker is accepted or not.\naccept_ratio::Float64: The acceptance ratio of the demon walk (excluding additional demon walks).\nat_final::AtomWalker: The final walker after the demon walks.\ndemon_energies::Array{Float64}: The energies of the demon at each step.\ntemp_estimate::Float64: The estimated temperature of the system.\n\n\n\n\n\n","category":"method"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.MC_random_walk!-Tuple{Int64, AtomWalker, LJParameters, Float64, Unitful.Quantity{Float64, 𝐋^2 𝐌 𝐓^-2, Unitful.FreeUnits{(eV,), 𝐋^2 𝐌 𝐓^-2, nothing}}}","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.MC_random_walk!","text":"MC_random_walk!(n_steps::Int, at::AtomWalker, lj::LJParameters, step_size::Float64, emax::typeof(0.0u\"eV\"))\n\nPerform a Monte Carlo random walk on the atomic/molecular system.\n\nArguments\n\nn_steps::Int: The number of Monte Carlo steps to perform.\nat::AtomWalker: The walker to perform the random walk on.\nlj::LJParameters: The Lennard-Jones potential parameters.\nstep_size::Float64: The maximum distance an atom can move in any direction.\nemax::typeof(0.0u\"eV\"): The maximum energy allowed for accepting a move.\n\nReturns\n\naccept_this_walker::Bool: Whether the walker is accepted or not.\naccept_rate::Float64: The acceptance rate of the random walk.\nat::AtomWalker: The updated walker.\n\n\n\n\n\n","category":"method"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.additional_demon_walk!-Tuple{Unitful.Quantity{Float64, 𝐋^2 𝐌 𝐓^-2, Unitful.FreeUnits{(eV,), 𝐋^2 𝐌 𝐓^-2, nothing}}, AtomWalker, LJParameters, Float64}","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.additional_demon_walk!","text":"additional_demon_walk!(e_demon::typeof(0.0u\"eV\"), at::AtomWalker, lj::LJParameters, step_size::Float64;\n                      e_demon_tolerance=1e-9u\"eV\", max_add_steps::Int=1_000_000)\n\nPerforms additional demon walk steps until the demon energy e_demon is below the tolerance e_demon_tolerance or the maximum number of additional steps max_add_steps is reached.\n\nArguments\n\ne_demon::typeof(0.0u\"eV\"): The initial demon energy.\nat::AtomWalker: The walker that the demon walk is performed on.\nlj::LJParameters: The LJ parameters.\nstep_size::Float64: The step size for the demon walk.\ne_demon_tolerance=1e-9u\"eV\": The tolerance for the demon energy.\nmax_add_steps::Int=1_000_000: The maximum number of additional steps.\n\nReturns\n\naccept_this_walker::Bool: Whether the walker is accepted after the additional demon walk.\naccept_rate::Float64: The acceptance rate of the additional demon walk.\nat::AtomWalker: The updated walker.\ne_demon::typeof(0.0u\"eV\"): The final demon energy.\n\n\n\n\n\n","category":"method"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.mean_sq_displacement-Tuple{AtomWalker, AtomWalker}","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.mean_sq_displacement","text":"mean_sq_displacement(at::AtomWalker, at_orig::AtomWalker)\n\nCalculate the mean squared displacement before and after random walk(s).\n\nArguments\n\nat::AtomWalker: The current AtomWalker after the random walk.\nat_orig::AtomWalker: The original AtomWalker before the random walk.\n\nReturns\n\ndistsq::typeof(0.0u\"Å\"^2): The mean squared displacement of all free particles.\n\n\n\n\n\n","category":"method"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.periodic_boundary_wrap!-Union{Tuple{T}, Tuple{StaticArraysCore.SVector{3, T}, AtomsBase.AbstractSystem}} where T","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.periodic_boundary_wrap!","text":"periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T\n\nWrap the position vector pos according to the periodic boundary conditions of the system. If the boundary condition is Periodic(), the position is wrapped using the modulo operator. If the boundary condition is DirichletZero(), the position is wrapped by reflecting the position vector across the boundary.\n\nArguments\n\npos::SVector{3,T}: The position vector to be wrapped.\nsystem::AbstractSystem: The system containing the periodic boundary conditions.\n\nReturns\n\nThe wrapped position vector.\n\n\n\n\n\n","category":"method"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.single_atom_demon_walk!-Tuple{AtomWalker, LJParameters, Float64}","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.single_atom_demon_walk!","text":"single_atom_demon_walk!(at::AtomWalker, lj::LJParameters, step_size::Float64;\n                        e_demon=0.0u\"eV\",\n                        demon_energy_threshold=Inf*u\"eV\",\n                        demon_gain_threshold=Inf*u\"eV\")\n\nPerform a single atom demon walk.\n\nArguments:\n\nat::AtomWalker: The walker to perform the demon walk on.\nlj::LJParameters: The LJ parameters.\nstep_size::Float64: The step size for the random walk.\n\nKeyword Arguments:\n\ne_demon=0.0u\"eV\": The energy of the demon.\ndemon_energy_threshold=Inf*u\"eV\": The energy threshold for the demon.\ndemon_gain_threshold=Inf*u\"eV\": The energy gain threshold for the demon during each move.\n\nReturns:\n\naccept::Bool: Whether the move is accepted or rejected.\nat::AtomWalker: The updated atom walker object.\ne_demon::Float64: The updated energy of the demon.\n\n\n\n\n\n","category":"method"},{"location":"AtomsMCMoves/#FreeBird.AtomsMCMoves.single_atom_random_walk!-Union{Tuple{T}, Tuple{StaticArraysCore.SVector{3, T}, Float64}} where T","page":"AtomsMCMoves","title":"FreeBird.AtomsMCMoves.single_atom_random_walk!","text":"single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64) where T\n\nPerform a single atom random walk by updating the position pos in each direction by a random amount. The step_size determines the maximum distance the atom can move in any direction.\n\nArguments\n\npos::SVector{3,T}: The current position of the atom as a 3D vector.\nstep_size::Float64: The maximum distance the atom can move in any direction.\n\nReturns\n\npos: The updated position of the atom.\n\n\n\n\n\n","category":"method"},{"location":"EnergyEval/#EnergyEval","page":"EnergyEval","title":"EnergyEval","text":"","category":"section"},{"location":"EnergyEval/#Functions","page":"EnergyEval","title":"Functions","text":"","category":"section"},{"location":"EnergyEval/","page":"EnergyEval","title":"EnergyEval","text":"Modules = [EnergyEval]","category":"page"},{"location":"EnergyEval/#FreeBird.EnergyEval","page":"EnergyEval","title":"FreeBird.EnergyEval","text":"Module for evaluating energy-related quantities for a system.\n\n\n\n\n\n","category":"module"},{"location":"EnergyEval/#FreeBird.EnergyEval.free_free_energy-Tuple{AtomsBase.AbstractSystem, LJParameters}","page":"EnergyEval","title":"FreeBird.EnergyEval.free_free_energy","text":"free_free_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)\n\nCalculate the energy from interactions between free particles using the Lennard-Jones potential. The energy is calculated by summing the pairwise interactions between the free particles.\n\nArguments\n\nat::AbstractSystem: The system for which the energy is calculated.\nlj::LJParameters: The Lennard-Jones parameters.\nfrozen::Int64: The number of frozen particles in the system.\n\nReturns\n\nfree_free_energy: The energy from interactions between free particles.\n\n\n\n\n\n","category":"method"},{"location":"EnergyEval/#FreeBird.EnergyEval.free_frozen_energy-Tuple{AtomsBase.AbstractSystem, LJParameters, Int64}","page":"EnergyEval","title":"FreeBird.EnergyEval.free_frozen_energy","text":"free_frozen_energy(at::AbstractSystem, lj::LJParameters, frozen::Int64)\n\nCompute the energy from interactions between free and frozen particles using the Lennard-Jones potential. The energy is calculated by summing the pairwise interactions between the free and frozen particles.\n\nArguments\n\nat::AbstractSystem: The system containing the particles.\nlj::LJParameters: The Lennard-Jones parameters.\nfrozen::Int64: The number of frozen particles.\n\nReturns\n\nfree_frozen_energy: The energy from interactions between free and frozen particles.\n\n\n\n\n\n","category":"method"},{"location":"EnergyEval/#FreeBird.EnergyEval.frozen_energy-Tuple{AtomsBase.AbstractSystem, LJParameters, Int64}","page":"EnergyEval","title":"FreeBird.EnergyEval.frozen_energy","text":"frozen_energy(at::AbstractSystem, lj::LJParameters, frozen::Int64)\n\nCompute the energy of the frozen particles in the system using the Lennard-Jones potential. The energy is calculated by summing the pairwise interactions between the frozen particles. Since the frozen particles do not move, the energy is typically only calculated once for a given system.\n\nArguments\n\nat::AbstractSystem: The system containing the particles.\nlj::LJParameters: The Lennard-Jones parameters.\nfrozen::Int64: The number of frozen particles.\n\nReturns\n\nfrozen_energy: The energy of the frozen particles in the system.\n\n\n\n\n\n","category":"method"},{"location":"EnergyEval/#FreeBird.EnergyEval.interaction_energy-Tuple{AtomsBase.AbstractSystem, LJParameters}","page":"EnergyEval","title":"FreeBird.EnergyEval.interaction_energy","text":"interaction_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)\n\nCompute the energy from interactions between particles using the Lennard-Jones potential. The energy is calculated by combining the energies from interactions between free particles  and between free and frozen particles. The energy from interactions between frozen particles is not included, as it is typically treated as a constant energy offset. If the number of frozen particles is zero, the energy calculated is equivalent to the total energy of the system.\n\nArguments\n\nat::AbstractSystem: The system containing the particles.\nlj::LJParameters: The Lennard-Jones parameters.\nfrozen::Int64: The number of frozen particles.\n\nReturns\n\ninteraction_energy: The energy from interactions between particles.\n\n\n\n\n\n","category":"method"},{"location":"EnergyEval/#FreeBird.EnergyEval.pbc_dist-Union{Tuple{T}, Tuple{Union{Vector{T}, StaticArraysCore.SVector{T}}, Union{Vector{T}, StaticArraysCore.SVector{T}}, AtomsBase.AbstractSystem}} where T","page":"EnergyEval","title":"FreeBird.EnergyEval.pbc_dist","text":"pbc_dist(pos1, pos2, at)\n\nCompute the distance between two positions considering periodic boundary conditions.\n\nArguments\n\npos1::Union{SVector{T},Vector{T}}: The first position.\npos2::Union{SVector{T},Vector{T}}: The second position.\nat::AbstractSystem: The abstract system containing boundary conditions and bounding box.\n\nReturns\n\nd::Float64: The distance between pos1 and pos2 considering periodic boundary conditions.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#EnergyEval","page":"FreeBirdIO","title":"EnergyEval","text":"","category":"section"},{"location":"FreeBirdIO/#Functions","page":"FreeBirdIO","title":"Functions","text":"","category":"section"},{"location":"FreeBirdIO/","page":"FreeBirdIO","title":"FreeBirdIO","text":"Modules = [FreeBirdIO]","category":"page"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO","text":"FreeBirdIO\n\nModule for input/output operations in the FreeBird package.\n\n\n\n\n\n","category":"module"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.DataSavingStrategy","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.DataSavingStrategy","text":"abstract type DataSavingStrategy\n\nAbstract type representing a strategy for saving data.\n\n\n\n\n\n","category":"type"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.SaveEveryN","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.SaveEveryN","text":"struct SaveEveryN <: DataSavingStrategy\n\nSaveEveryN is a concrete subtype of DataSavingStrategy that specifies saving data every N steps.\n\nFields\n\nfilename::String: The name of the file to save the data to.\nn::Int: The number of steps between each save.\n\n\n\n\n\n","category":"type"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.convert_system_to_walker-Tuple{AtomsBase.FlexibleSystem, Bool}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.convert_system_to_walker","text":"convert_system_to_walker(at::FlexibleSystem, resume::Bool)\n\nConverts a FlexibleSystem object to an AtomWalker object.\n\nArguments\n\nat::FlexibleSystem: The FlexibleSystem object to convert.\nresume::Bool: Whether to resume from previous data.\n\nReturns\n\nAtomWalker: The converted AtomWalker object.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.convert_walker_to_system-Tuple{AtomWalker}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.convert_walker_to_system","text":"convert_walker_to_system(at::AtomWalker)\n\nConverts an AtomWalker object to an AbstractSystem object.\n\nArguments\n\nat::AtomWalker: The AtomWalker object to be converted.\n\nReturns\n\nAbstractSystem: The converted AbstractSystem object.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.generate_initial_configs-Tuple{Int64, Float64, Int64}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.generate_initial_configs","text":"generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)\n\nGenerate initial configurations for a given number of walkers.\n\nArguments\n\nnum_walkers::Int: The number of walkers.\nvolume_per_particle::Float64: The volume per particle.\nnum_particle::Int: The number of particles.\nparticle_type::Symbol=:H: The type of particle (default is :H).\n\nReturns\n\nAn array of initial configurations for each walker.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.generate_random_starting_config-Tuple{Float64, Int64}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.generate_random_starting_config","text":"generate_random_starting_config(volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)\n\nGenerate a random starting configuration for a system of particles.\n\nArguments\n\nvolume_per_particle::Float64: The volume per particle.\nnum_particle::Int: The number of particles.\nparticle_type::Symbol=:H: The type of particle (default is hydrogen).\n\nReturns\n\nFastSystem: A FastSystem object representing the generated system.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.read_configs-Tuple{String, Vector}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.read_configs","text":"read_configs(filename::String, pbc::Vector)\n\nReads atomic configurations from a file and applies periodic boundary conditions.\n\nArguments\n\nfilename::String: The name of the file containing the atomic configurations.\npbc::Vector: A vector specifying the periodic boundary conditions.\n\nReturns\n\nAn array of atomic configurations with periodic boundary conditions applied.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.read_configs-Tuple{String}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.read_configs","text":"read_configs(filename::String; pbc::String=\"TTT\")\n\nReads configurations from a file.\n\nArguments\n\nfilename::String: The name of the file to read configurations from.\npbc::String=\"TTT\": Periodic boundary conditions. A string of length 3, where each character represents whether the corresponding dimension has periodic boundary conditions ('T') or not ('F').\n\nReturns\n\nThe configurations read from the file.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.read_single_config-Tuple{String, Vector}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.read_single_config","text":"read_single_config(filename::String, pbc::Vector)\n\nReads a single configuration from the specified file and sets the periodic boundary conditions (PBC) for the atoms.\n\nArguments\n\nfilename::String: The name of the file to read the configuration from.\npbc::Vector: A vector specifying the periodic boundary conditions.\n\nReturns\n\nat::Atoms: The atoms with the PBC set.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.read_single_config-Tuple{String}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.read_single_config","text":"read_single_config(filename::String; pbc::String=\"TTT\")\n\nReads a single configuration from the specified file.\n\nArguments\n\nfilename::String: The name of the file to read from.\npbc::String=\"TTT\": The periodic boundary conditions. Default is \"TTT\".\n\nReturns\n\nThe configuration read from the file.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.read_single_walker-Tuple{String}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.read_single_walker","text":"read_single_walker(filename::String; pbc::String=\"TTT\", resume::Bool=true)\n\nReads a single walker from the specified file.\n\nArguments\n\nfilename::String: The path to the file containing the walker data.\npbc::String: (optional) The periodic boundary conditions. Default is \"TTT\".\nresume::Bool: (optional) Whether to resume reading from a previous checkpoint. Default is true.\n\nReturns\n\nThe walker object read from the file.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.read_walkers-Tuple{String}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.read_walkers","text":"read_walkers(filename::String; pbc::String=\"TTT\", resume::Bool=true)\n\nReads walker configurations from a file.\n\nArguments\n\nfilename::String: The name of the file to read the walker configurations from.\npbc::String: A string specifying the periodic boundary conditions. Default is \"TTT\".\nresume::Bool: A boolean indicating whether to resume reading from a previous checkpoint. Default is true.\n\nReturns\n\nAn array of walker objects.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.set_pbc-Tuple{ExtXYZ.Atoms, Vector}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.set_pbc","text":"set_pbc(at::Atoms, pbc::Vector)\n\nSet the periodic boundary conditions for a system of atoms.\n\nArguments\n\nat::Atoms: The system of atoms.\npbc::Vector: A vector of length 3 specifying the periodic boundary conditions for each dimension. Each element can be either true for periodic boundary conditions or false for Dirichlet zero boundary conditions.\n\nReturns\n\nFlexibleSystem: A flexible system with the specified boundary conditions.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.write_df-Tuple{String, DataFrames.DataFrame}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.write_df","text":"write_df(filename::String, df::DataFrame)\n\nWrite a DataFrame to a CSV file.\n\nArguments\n\nfilename::String: The name of the file to write to.\ndf::DataFrame: The DataFrame to write.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.write_df_every_n-Tuple{DataFrames.DataFrame, Int64, SaveEveryN}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.write_df_every_n","text":"write_df_every_n(df::DataFrame, step::Int, d_strategy::SaveEveryN)\n\nWrite the DataFrame df to a file specified by d_strategy.filename every d_strategy.n steps.\n\nArguments\n\ndf::DataFrame: The DataFrame to be written.\nstep::Int: The current step number.\nd_strategy::SaveEveryN: The save strategy specifying the filename and the step interval.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.write_single_walker-Tuple{String, AtomWalker}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.write_single_walker","text":"write_single_walker(filename::String, at::AtomWalker)\n\nWrite a single AtomWalker object to a file.\n\nArguments\n\nfilename::String: The name of the file to write to.\nat::AtomWalker: The AtomWalker object to write.\n\n\n\n\n\n","category":"method"},{"location":"FreeBirdIO/#FreeBird.FreeBirdIO.write_walkers-Tuple{String, Vector{AtomWalker}}","page":"FreeBirdIO","title":"FreeBird.FreeBirdIO.write_walkers","text":"write_walkers(filename::String, ats::Vector{AtomWalker})\n\nWrite a collection of AtomWalker objects to a file.\n\nArguments\n\nfilename::String: The name of the file to write the walkers to.\nats::Vector{AtomWalker}: The collection of AtomWalker objects to write.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#Sampling-Schemes","page":"SamplingSchemes","title":"Sampling Schemes","text":"","category":"section"},{"location":"SamplingSchemes/#Functions","page":"SamplingSchemes","title":"Functions","text":"","category":"section"},{"location":"SamplingSchemes/","page":"SamplingSchemes","title":"SamplingSchemes","text":"Modules = [SamplingSchemes]","category":"page"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.MCDemonWalk","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.MCDemonWalk","text":"MCDemonWalk <: MCRoutine\n\nA Monte Carlo routine for performing demon walks.\n\nFields\n\ne_demon_tolerance::typeof(0.0u\"eV\"): The tolerance for the energy difference in the demon walk.\ndemon_energy_threshold::typeof(0.0u\"eV\"): The energy threshold for the demon walk.\ndemon_gain_threshold::typeof(0.0u\"eV\"): The gain threshold for the demon during each walk.\nmax_add_steps::Int: The maximum number of steps to add in the demon walk if the demon energy e_demon is higher than the demon_energy_threshold.\n\n\n\n\n\n","category":"type"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.MCRandomWalk","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.MCRandomWalk","text":"struct MCRandomWalk <: MCRoutine\n\nA type representing a Monte Carlo random walk sampling scheme.\n\n\n\n\n\n","category":"type"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.MCRoutine","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.MCRoutine","text":"abstract type MCRoutine\n\nAn abstract type representing a Monte Carlo routine.\n\n\n\n\n\n","category":"type"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.MixedMCRoutine","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.MixedMCRoutine","text":"struct MixedMCRoutine <: MCRoutine\n\nA mutable struct representing a mixed Monte Carlo routine, where the main routine is used for the majority of the steps,      and the backup routine is used when the main routine fails to accept a move. Currently, it is intended to use MCRandomWalk      as the main routine and MCDemonWalk as the backup routine.\n\nFields\n\nmain_routine::MCRoutine: The main Monte Carlo routine.\nback_up_routine::MCRoutine: The backup Monte Carlo routine.\nns_params_main::NestedSamplingParameters: The nested sampling parameters for the main routine.\nns_params_back_up::NestedSamplingParameters: The nested sampling parameters for the backup routine.\n\n\n\n\n\n","category":"type"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.NestedSamplingParameters","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.NestedSamplingParameters","text":"mutable struct NestedSamplingParameters <: SamplingParameters\n\nThe NestedSamplingParameters struct represents the parameters used in the nested sampling scheme.\n\nFields\n\nmc_steps::Int64: The number of total Monte Carlo moves to perform.\nstep_size::Float64: The step size used in the sampling process.\nfail_count::Int64: The number of failed MC moves in a row. Used to terminate the sampling process if it exceeds a certain threshold.\n\n\n\n\n\n","category":"type"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.adjust_step_size-Tuple{NestedSamplingParameters, Float64}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.adjust_step_size","text":"adjust_step_size(ns_params::NestedSamplingParameters, rate::Float64)\n\nAdjusts the step size of the nested sampling algorithm based on the acceptance rate.      If the acceptance rate is greater than 0.75, the step size is increased by 1%.      If the acceptance rate is less than 0.25, the step size is decreased by 1%.\n\nArguments\n\nns_params::NestedSamplingParameters: The parameters of the nested sampling algorithm.\nrate::Float64: The acceptance rate of the algorithm.\n\nReturns\n\nns_params::NestedSamplingParameters: The updated parameters with adjusted step size.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.nested_sampling_loop!-Tuple{AtomWalkers, Int64, MixedMCRoutine, DataSavingStrategy}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.nested_sampling_loop!","text":"nested_sampling_loop!(liveset::AtomWalkers, n_steps::Int64, mc_routine::MixedMCRoutine, save_strategy::DataSavingStrategy)\n\nPerform a nested sampling loop for a given number of steps.\n\nArguments\n\nliveset::AtomWalkers: The initial set of walkers.\nn_steps::Int64: The number of steps to perform.\nmc_routine::MixedMCRoutine: The mixed Monte Carlo routine to use.\nsave_strategy::DataSavingStrategy: The strategy for saving data.\n\nReturns\n\ndf::DataFrame: The data frame containing the iteration number and maximum energy for each step.\nliveset::AtomWalkers: The updated set of walkers.\nmc_routine.ns_params_main: The updated nested sampling parameters for the main routine.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.nested_sampling_loop!-Tuple{AtomWalkers, NestedSamplingParameters, Int64, MCRoutine}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.nested_sampling_loop!","text":"nested_sampling_loop!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, n_steps::Int64, mc_routine::MCRoutine; args...)\n\nPerform a nested sampling loop for a given number of steps.\n\nArguments\n\nliveset::AtomWalkers: The initial set of walkers.\nns_params::NestedSamplingParameters: The parameters for nested sampling.\nn_steps::Int64: The number of steps to perform.\nmc_routine::MCRoutine: The Monte Carlo routine to use.\n\nKeyword Arguments\n\nargs...: Additional arguments.\n\nReturns\n\ndf: A DataFrame containing the iteration number and maximum energy for each step.\nliveset: The updated set of walkers.\nns_params: The updated nested sampling parameters.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.nested_sampling_step!-Tuple{AtomWalkers, NestedSamplingParameters, MCDemonWalk}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.nested_sampling_step!","text":"nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCDemonWalk)\n\nPerform a single step of the nested sampling algorithm using the Monte Carlo demon walk routine.\n\nArguments\n\nliveset::AtomWalkers: The set of atom walkers representing the current state of the system.\nns_params::NestedSamplingParameters: The parameters for the nested sampling algorithm.\nmc_routine::MCDemonWalk: The parameters for the Monte Carlo demon walk.\n\nReturns\n\niter: The iteration number after the step.\nemax: The maximum energy recorded during the step.\nliveset: The updated set of atom walkers after the step.\nns_params: The updated nested sampling parameters after the step.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.nested_sampling_step!-Tuple{AtomWalkers, NestedSamplingParameters, MCRandomWalk}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.nested_sampling_step!","text":"nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRandomWalk)\n\nPerform a single step of the nested sampling algorithm using the Monte Carlo random walk routine.\n\nArguments\n\nliveset::AtomWalkers: The set of atom walkers.\nns_params::NestedSamplingParameters: The parameters for nested sampling.\nmc_routine::MCRandomWalk: The Monte Carlo random walk routine.\n\nReturns\n\niter: The iteration number after the step.\nemax: The highest energy recorded during the step.\nliveset: The updated set of atom walkers.\nns_params: The updated nested sampling parameters.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.sort_by_energy!-Tuple{AtomWalkers}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.sort_by_energy!","text":"sort_by_energy!(liveset::LJAtomWalkers)\n\nSorts the walkers in the liveset by their energy in descending order.\n\nArguments\n\nliveset::LJAtomWalkers: The liveset of walkers to be sorted.\n\nReturns\n\nliveset::LJAtomWalkers: The sorted liveset.\n\n\n\n\n\n","category":"method"},{"location":"SamplingSchemes/#FreeBird.SamplingSchemes.update_iter!-Tuple{AtomWalkers}","page":"SamplingSchemes","title":"FreeBird.SamplingSchemes.update_iter!","text":"update_iter!(liveset::AtomWalkers)\n\nUpdate the iteration count for each walker in the liveset.\n\nArguments\n\nliveset::AtomWalkers: The set of walkers to update.\n\n\n\n\n\n","category":"method"},{"location":"Potentials/#Potentials","page":"Potentials","title":"Potentials","text":"","category":"section"},{"location":"Potentials/#Functions","page":"Potentials","title":"Functions","text":"","category":"section"},{"location":"Potentials/","page":"Potentials","title":"Potentials","text":"Modules = [Potentials]","category":"page"},{"location":"Potentials/#FreeBird.Potentials","page":"Potentials","title":"FreeBird.Potentials","text":"Module for defining and implementing potentials.\n\n\n\n\n\n","category":"module"},{"location":"Potentials/#FreeBird.Potentials.LJParameters","page":"Potentials","title":"FreeBird.Potentials.LJParameters","text":"struct LJParameters\n\nThe LJParameters struct represents the parameters for the Lennard-Jones potential.\n\nFields\n\nepsilon::typeof(1.0u\"eV\"): The energy scale of the potential.\nsigma::typeof(1.0u\"Å\"): The length scale of the potential.\ncutoff::Float64: The cutoff distance for the potential, in units of sigma.\nshift::typeof(0.0u\"eV\"): The energy shift applied to the potential, calculated at the cutoff distance.\n\n\n\n\n\n","category":"type"},{"location":"Potentials/#FreeBird.Potentials.LJParameters-Tuple{}","page":"Potentials","title":"FreeBird.Potentials.LJParameters","text":"LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)\n\nA constructor for the LJParameters struct with default values for the  Lennard-Jones potential with no cutoff or shift. The shift parameter can be  specified as a boolean, if true, the shift energy is calculated automatically  at the cutoff distance; or as a typeof(0.0u\"eV\"), in which case the value is used directly.\n\nExample\n\njulia> lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)\nLJParameters(0.1 eV, 2.5 Å, 3.5, 0.0 eV)\n\njulia> lj = LJParameters(sigma=2.5)\nLJParameters(1.0 eV, 2.5 Å, Inf, 0.0 eV)\n\njulia> lj = LJParameters(cutoff=3.5,shift=5.0)\nLJParameters(1.0 eV, 1.0 Å, 3.5, 5.0 eV)\n\njulia> lj = LJParameters(cutoff=3.5,shift=true)\nLJParameters(1.0 eV, 1.0 Å, 3.5, -0.0021747803916549904 eV)\n\njulia> lj = LJParameters(cutoff=3.5,shift=false)\nLJParameters(1.0 eV, 1.0 Å, 3.5, 0.0 eV)\n\n\n\n\n\n\n","category":"method"},{"location":"Potentials/#FreeBird.Potentials.lj_energy-Tuple{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}, LJParameters}","page":"Potentials","title":"FreeBird.Potentials.lj_energy","text":"lj_energy(r::typeof(1.0u\"Å\"), lj::LJParameters)\n\nCompute the Lennard-Jones energy between two particles at a given distance.\n\nArguments\n\nr::typeof(1.0u\"Å\"): The distance between the particles.\nlj::LJParameters: The Lennard-Jones parameters.\n\nReturns\n\n0.0u\"eV\" if the distance is greater than the cutoff distance.\nThe Lennard-Jones energy minus the shift otherwise.\n\n\n\n\n\n","category":"method"},{"location":"Potentials/#FreeBird.Potentials.lj_energy-Tuple{Unitful.Quantity{Float64, 𝐋^2 𝐌 𝐓^-2, Unitful.FreeUnits{(eV,), 𝐋^2 𝐌 𝐓^-2, nothing}}, Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}, Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}","page":"Potentials","title":"FreeBird.Potentials.lj_energy","text":"lj_energy(epsilon::typeof(1.0u\"eV\"), sigma::typeof(1.0u\"Å\"), r::typeof(1.0u\"Å\"))\n\nCompute the Lennard-Jones potential energy between two particles.\n\nThe Lennard-Jones potential energy is given by the equation:\n\nV(r_ij) = 4varepsilon_ij leftleft(fracsigma_ijr_ijright)^12 - left(fracsigma_ijr_ijright)^6right\n\nwhere epsilon is the energy scale, sigma is the distance scale, and r is the distance between the particles.\n\nArguments\n\nepsilon::typeof(1.0u\"eV\"): The energy scale of the potential.\nsigma::typeof(1.0u\"Å\"): The distance scale of the potential.\nr::typeof(1.0u\"Å\"): The distance between the particles.\n\nReturns\n\nThe Lennard-Jones potential energy between the particles.\n\n\n\n\n\n","category":"method"},{"location":"AbstractWalkers/#AbstractWalkers","page":"AbstractWalkers","title":"AbstractWalkers","text":"","category":"section"},{"location":"AbstractWalkers/#Functions","page":"AbstractWalkers","title":"Functions","text":"","category":"section"},{"location":"AbstractWalkers/","page":"AbstractWalkers","title":"AbstractWalkers","text":"Modules = [AbstractWalkers]","category":"page"},{"location":"AbstractWalkers/#FreeBird.AbstractWalkers","page":"AbstractWalkers","title":"FreeBird.AbstractWalkers","text":"Module containing abstract definitions for walkers.\n\n\n\n\n\n","category":"module"},{"location":"AbstractWalkers/#FreeBird.AbstractWalkers.AtomWalker","page":"AbstractWalkers","title":"FreeBird.AbstractWalkers.AtomWalker","text":"mutable struct AtomWalker\n\nThe AtomWalker struct represents a walker composed of atoms/molecules.\n\nFields\n\nconfiguration::FastSystem: The configuration of the walker.\nenergy::typeof(0.0u\"eV\"): The energy of the walker.\niter::Int64: The current iteration number of the walker.\nnum_frozen_part::Int64: The number of frozen particles in the walker.\nenergy_frozen_part::typeof(0.0u\"eV\"): The energy of the frozen particles in the walker, serves as a constant energy offset   to the interacting part of the system.\n\nConstructor\n\nAtomWalker(configuration::FastSystem; energy=0.0u\"eV\", iter=0, num_frozen_part=0, energy_frozen_part=0.0u\"eV\")\n\nCreate a new AtomWalker with the given configuration and optional energy, iteration number, number of frozen particles, and energy of the frozen particles.\n\n\n\n\n\n","category":"type"},{"location":"AbstractWalkers/#FreeBird.AbstractWalkers.LJAtomWalkers","page":"AbstractWalkers","title":"FreeBird.AbstractWalkers.LJAtomWalkers","text":"struct LJAtomWalkers <: AtomWalkers\n\nThe LJAtomWalkers struct contains a list of AtomWalker objects and the Lennard-Jones potential parameters      that defines the interactions between the particles in the walkers.\n\nFields\n\nwalkers::Vector{AtomWalker}: The list of AtomWalkers.\nlj_potential::LJParameters: The Lennard-Jones potential parameters.\n\nConstructors\n\nLJAtomWalkers(walkers::Vector{AtomWalker}, lj_potential::LJParameters): Constructs a new LJAtomWalkers    object with the given walkers and Lennard-Jones potential parameters. The energies of the walkers are automatically   assigned using the Lennard-Jones potential parameters.\n\n\n\n\n\n","category":"type"},{"location":"AbstractWalkers/#FreeBird.AbstractWalkers.LJAtomWalkers-Tuple{Vector{AtomWalker}, LJParameters, Int64}","page":"AbstractWalkers","title":"FreeBird.AbstractWalkers.LJAtomWalkers","text":"LJAtomWalkers(ats::Vector{AtomWalker}, lj::LJParameters, num_frozen_part::Int64)\n\nCreate a new LJAtomWalkers object with the given number of frozen particles. This function is a convenient  wrapper around the LJAtomWalkers constructor when the number of frozen particles is undefined or modified.\n\nArguments\n\nats::Vector{AtomWalker}: The list of AtomWalker objects.\nlj::LJParameters: The Lennard-Jones potential parameters.\nnum_frozen_part::Int64: The number of frozen particles.\n\nReturns\n\nLJAtomWalkers: The LJAtomWalkers object with the updated number of frozen particles and energies assigned to the walkers.\n\n\n\n\n\n","category":"method"},{"location":"AbstractWalkers/#FreeBird.AbstractWalkers.assign_lj_energies!-Tuple{Vector{AtomWalker}, LJParameters}","page":"AbstractWalkers","title":"FreeBird.AbstractWalkers.assign_lj_energies!","text":"assign_lj_energies!(walkers::Vector{AtomWalker}, lj::LJParameters)\n\nAssigns LJ energies to a list of AtomWalkers.\n\nArguments\n\nwalkers::Vector{AtomWalker}: A list of AtomWalkers to assign energies to.\nlj::LJParameters: The LJ parameters used for energy calculation.\n\nReturns\n\nwalkers::Vector{AtomWalker}: The updated list of AtomWalkers with assigned energies.\n\n\n\n\n\n","category":"method"},{"location":"AbstractWalkers/#FreeBird.AbstractWalkers.update_walker!-Tuple{AtomWalker, Symbol, Any}","page":"AbstractWalkers","title":"FreeBird.AbstractWalkers.update_walker!","text":"update_walker!(walker::AtomWalker, key::Symbol, value)\n\nUpdate the properties of an AtomWalker object.\n\nA convenient function that updates the value of a specific property of an AtomWalker object.\n\nArguments\n\nwalker::AtomWalker: The AtomWalker object to be updated.\nkey::Symbol: The key of the property to be updated.\nvalue: The new value of the property.\n\nReturns\n\nwalker::AtomWalker: The updated AtomWalker object.\n\nExample\n\nupdate_walker!(walker, :energy, 10.0u\"eV\")\nupdate_walker!(walker, :iter, 1)\n\n\n\n\n\n","category":"method"},{"location":"#FreeBird.jl","page":"FreeBird.jl","title":"FreeBird.jl","text":"","category":"section"},{"location":"","page":"FreeBird.jl","title":"FreeBird.jl","text":"Documentation for FreeBird.jl","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"","category":"page"}]
}