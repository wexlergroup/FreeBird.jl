# # Quick start guide to FreeBird.jl

# This is a quick start guide to using the FreeBird.jl package.
# It covers the basic functionalities of the package, such as 
# generating atomistic and lattice walkers, defining a potential energy function or
# Hamiltonian, and running a sampling simulation.
# For more detailed information, please refer to the documentation of the package.
# You can find the runnable version of this script in the `scripts` directory of the package.

# ## Atomistic walkers and nested sampling

# First, let's load the FreeBird.jl package:
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
# Here, `ls` is a [`LJAtomWalkers`](@ref) type, and has the `walkers` and `lj` fields attached to it.

# Now, time to set up a simulation. We will be using nested sampling, a Bayesian-inference inspired method, as an example here.
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
save = SaveEveryN(n_traj=10, n_snap=20_000)
# The [`SaveEveryN`](@ref) type is a struct that holds the parameters of the saving routine. 

# Now, we are ready to run the nested sampling simulation:
energies, liveset, _ = nested_sampling_loop!(ls, ns_params, 20_000, mc, save)
# The results of the simulation are stored in the `energies` and `liveset` variables.
# The `energies` variable is a `DataFrame` that contains the energies of the walkers at each iteration.
# The `liveset` variable is the final liveset after the simulation.
# Let's see how the walkers look like after the simulation:
liveset.walkers[1].configuration
# They should be in a more ordered state, in this case, a cluster, than the initial gaseous state.

# ### Calculating heat capacity with `AnalysisTools` module

# The [`AnalysisTools`](@ref) module provides functions to calculate the heat capacity of the system.
# First, we calculate the `ω` factors, which account for the fractions of phase-space volume sampled during
# each nested sampling iteration, defined as:
# ```math
# \omega_i = \frac{1}{N+1} \left(\frac{N}{N+1}\right)^i
# ````
# where $N$ is the number of walkers and $i$ is the iteration number.
ωi = ωᵢ(energies.iter, 120)
# Let's shift the energies to be greater than or equal to zero, making the calculation of the heat capacity more stable.
Ei = energies.emax .- minimum(energies.emax)
# Specify the temperatures that we are interested in, in units of Kelvin.
Ts = collect(1:0.1:1000)
# Define the Boltzmann constant in units of eV/K.
kb = 8.617333262e-5 # eV/K
# Calculate the inverse temperatures
β = 1 ./(kb.*Ts)
# Define the degrees of freedom, which is 3×6 for the 6-particle system.
dof = 18
# Calculate the heat capacities as a function of temperature using the `cv` function,
# ```math
# C_V(\beta) = \frac{\mathrm{dof} \cdot k_B}{2} + k_B \beta^2 \left(\frac{\sum_i \omega_i E_i^2 \exp(-E_i \beta)}{Z(\beta)} - U(\beta)^2\right)
# ```
cvs = cv(energies, β, dof, 120)

# Let's plot the heat capacity as a function of temperature
using Plots
plot(Ts, cvs./kb, xlabel="Temperature (K)", ylabel="Heat Capacity (\$k_B\$)", label="LJ\$_6\$")

# The plot should show the heat capacity as a function of temperature for the 6-particle Lennard-Jones system, with a main peak around 400 K, representing the phase transition, and some fluctuations at low temperatures, and tailing off to zero at high temperatures.

# That's it! You have successfully run a nested sampling simulation using the FreeBird.jl package.

# ## Lattice walkers and exact enumeration

# Another feature of FreeBird.jl is the ability to work with lattice systems.
# The lattice systems are defined by the [`MLattice`](@ref) which is a parametrized type.
#=====================
```julia
MLattice{C,G}(
    lattice_vectors::Matrix{Float64},
    basis::Vector{Tuple{Float64, Float64, Float64}},
    supercell_dimensions::Tuple{Int64, Int64, Int64},
    periodicity::Tuple{Bool, Bool, Bool},
    components::Vector{Vector{Bool}},
    adsorptions::Vector{Bool},
    cutoff_radii::Vector{Float64},
) where {C,G}
```
=====================#
# The `C` parameter is the number of components in the system, and the `G` parameter defines the geometry of the lattice.

# Now, let's create a simple square lattice system with single component:
ml = MLattice{1,SquareLattice}(components=[[1,2,3,4]])

#==
When you run the above code, the outer constructor of `MLattice` will be called.
Many of the arguments are optional and have default values.
The `components` argument is a vector of vectors that defines the components of the system.
The `components=[[1,2,3,4]]` argument specifies that the system has a single component,
and the first four sites are occupied.
```julia
MLattice{C,SquareLattice}(; lattice_constant::Float64=1.0,
    basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0)],
    supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 4, 1),
    periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
    cutoff_radii::Vector{Float64}=[1.1, 1.5],
    components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
    adsorptions::Union{Vector{Int},Symbol}=:full)
```
==#

# You may notice that the above code returns a `SLattice` type.
# The `SLattice` type is simply an alias for the `MLattice{1,G}`,
# where `G` is the geometry of the lattice and the number of components is fixed to 1.
# You can also directly call the `SLattice`, it will give the same result:
sl = SLattice{SquareLattice}(components=[[1,2,3,4]])

# Now, let's define a Hamiltonian for the lattice system:
ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
# The [`GenericLatticeHamiltonian`](@ref) type is a struct that holds the parameters of the Hamiltonian.
# The first argument is the on-site energy, and the second argument is the list of n-th nearest-neighbors energy.
# The third argument is the unit of the energy.

# To run exact enumeration, we only need a initial walker/lattice configuration, and
# the Hamiltonian. Let's run the exact enumeration:
df, ls = exact_enumeration(sl, ham)

# The results of the exact enumeration are stored in the `df` and `ls` variables.
# The `df` variable is a `DataFrame` that contains the list of energies, as well as the configurations.
# The `ls` variable is the final liveset that contains all possible configurations of the lattice system.
# Let's see how the first configuration looks like:
ls.walkers[1].configuration
# It's the initial configuration of the lattice system.
# Let's see how the last configuration looks like:
ls.walkers[end].configuration

# Be warned that the exact enumeration can be computationally expensive for large systems.

# ### Calculating heat capacity

# Since we enumerated all possible configurations of the lattice system, we can calculate the partition function, then heat capacity directly.

# Let's calculate the heat capacity for the lattice system:
# Define the temperatures that we are interested in, in units of Kelvin.
Ts = collect(1:0.1:500)
# Convert them to inverse temperatures
βs = 1 ./(kb.*Ts)
# Extract the energies from the DataFrame, keeping the values only
es = [e.val for e in df.energy]

# Since this is not a nested sampling run, each configuration carries the same weight:
ω_1 = ones(length(df.energy))
# And for a lattice, the degrees of freedom is 0:
dof = 0
# Now we can use a scaler version of the [`cv`](@ref) function to calculate the heat capacity:
cvs = [cv(β, ω_1, es, dof) for β in βs]

# Let's plot the heat capacity as a function of temperature
plot(Ts, cvs./kb, xlabel="Temperature (K)", ylabel="Heat Capacity (\$k_B\$)", label="Square Lattice")
# You should expect to see a single peak in the heat capacity curve around 40 K, and tailing off to zero at high temperatures.

# That's it! You have successfully run an exact enumeration simulation using the FreeBird.jl package.