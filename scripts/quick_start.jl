# This is a quick start guide to using the FreeBird.jl package.

# Load the FreeBird.jl package:
using FreeBird

# Now, let's create a few configurations of a simple atomistic system with six particles in a 3D box.
single_config = generate_initial_configs(1, 562.5, 6; particle_type=:H)

# The function above has generated a single configuration, with 562.5 Å$^3$ volume per particle, and 6 particles of type H.
# Note that the `particle_type` keyword argument can be used to specify the type of particle, i.e., chemical element. By default, the type is set to `:H`.
# Use `?generate_initial_configs` in the REPL to see the documentation of the function. Or see [`generate_initial_configs`](@ref).

# Let's inspect the generated configuration:
single_config[1]

# It's of a `FastSystem` type from `AtomsBase`. The dimensions of the box are 15 Å x 15 Å x 15 Å, following the volume per particle specified.
# The positions of the particles are randomly generated within the box.

# Now, let's generate a few more configurations:
configs = generate_initial_configs(120, 562.5, 6)

# The function above has generated 120 configurations, again, with 562.5 Å$^3$ volume per particle, and 6 particles of type H.
# These configurations will be served as the initial walkers for a sampling run, but first, we need to warp them into the [`AtomWalker`](@ref) type defined in FreeBird.jl.
walkers = AtomWalker.(generate_initial_configs(120, 562.5, 6))

# Let's inquire the type of the `walkers` variable:
walkers |> typeof

# The `walkers` variable is of a `Vector{AtomWalker{1}}` type, which is a vector of `AtomWalker{1}` objects.
# The `AtomWalker{1}` type is a parametrized type, where the parameter is the number of components in the system.
# In this case, the system has only one component, consisting of 6 particles of type H.

# To define how these particles interact with each other, we need to create a potential energy function. Let's use the Lennard-Jones potential:
lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0)

# The [`LJParameters`](@ref) type is a struct that holds the parameters of the Lennard-Jones potential.
# The `epsilon` and `sigma` fields are the energy and length scales of the potential, respectively.
# The `cutoff` field is the distance at which the potential is truncated. 
# An energy shift is applied to the potential to ensure continuity at the cutoff distance, automatically in this case.
# See the documentation of the `LJParameters` type for more information and examples.

# We now can create a so-called *liveset* that will be used to store the walkers during the simulation. 
# The `lj` potential will be used to attached and used to calculate the potential energy of the walkers.
ls = LJAtomWalkers(walkers, lj)
# Here, `ls` is a [`LJAtomWalkers`](@ref) type, and has the `walkers` and `lj` files attached to it.

# Now, time to set up a simulation. We will be using nested sampling, a Bayesian inference method, as an example here.
# First, we need to define the nested sampling parameters:
ns_params = NestedSamplingParameters(200, 0.1, 0.01, 1e-5, 1.0, 0, 200)

# The [`NestedSamplingParameters`](@ref) type is a struct that holds the parameters of the nested sampling algorithm.
# The fields are as follows:
# -  `mc_steps::Int64`: The number of total Monte Carlo moves to perform.
# -  `initial_step_size`::Float64: The initial step size, which is the fallback step size if MC routine fails to accept a move.
# -  `step_size::Float64`: The on-the-fly step size used in the sampling process.
# -  `step_size_lo::Float64`: The lower bound of the step size.
# -  `step_size_up::Float64`: The upper bound of the step size.
# -  `fail_count::Int64`: The number of failed MC moves in a row.
# -  `allowed_fail_count::Int64`: The maximum number of failed MC moves allowed before resetting the step size.

# Speaking of the Monte Carlo moves, we need to define that too:
mc = MCRandomWalkClone()
# The [`MCRandomWalkClone`](@ref) type is a type of Monte Carlo move that indicates a new walker is created by cloning an existing walker and then decorrelate the positions of the particles.

# We also need to specify how we want to save the data and the output:
save = SaveEveryN(n_traj=10, n_snap=20_00)
# The [`SaveEveryN`](@ref) type is a struct that holds the parameters of the saving routine. 

# Now, we are ready to run the nested sampling simulation:
energies, liveset, _ = nested_sampling_loop!(ls, ns_params, 20_000, mc, save)
# The results of the simulation are stored in the `energies` and `liveset` variables.
# The `energies` variable is a `DataFrame` that contains the energies of the walkers at each iteration.
# The `liveset` variable is the final liveset after the simulation.
# Let's see how the walkers look like after the simulation:
liveset.walkers[1].configuration
# They should be in a more ordered state, in this case, a cluster, than the initial gaseous state.

# That's it! You have successfully run a nested sampling simulation using the FreeBird.jl package.

# For more information, please refer to the documentation of the FreeBird.jl package.

# Please see other scripts in the `scripts` directory for more examples and use cases.
# More tutorials and examples will be added in the future. Stay tuned!