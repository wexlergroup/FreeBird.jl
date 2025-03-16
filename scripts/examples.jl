# # FreeBird.jl Examples

# Here we provide a few examples of how to use FreeBird.jl. 
# We focus on how to set up a sampling run with various sampling schemes.
# You may need to fine-tune the parameters for your specific system.

# ## Atomistic system

# Before we start, we need to load the FreeBird.jl package.

using FreeBird

# ### Nested Sampling

# First, let's generate some initial configurations
configs = generate_initial_configs(120, 562.5, 6)

# Walkers are the objects that will be used in the simulation.
# They contain the information about the system's configuration.
# Let's warp the configurations into walkers.
walkers = AtomWalker.(generate_initial_configs(120, 562.5, 6))
# Let's define the Lennard-Jones potential
lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0)
# We then create a `LJAtomWalkers` object that contains the walkers and the Lennard-Jones potential.
ls = LJAtomWalkers(walkers, lj)

# We then set up a `NestedSamplingParameters` object with the desired parameters. See the documentation for more information.
ns_params = NestedSamplingParameters(mc_steps=200, step_size=0.1, random_seed=1234*rand(Int))
# We then set up how we want to perform the Monte Carlo moves. `MCRandomWalkClone` is a routine that clones an existing walker and performs random walks.
mc = MCRandomWalkClone()
# We also set up a `SaveEveryN` object for saving the output.
save = SaveEveryN(n_traj=10, n_snap=20_000, df_filename="output_df_lj6.csv", n_info=10)
# Finally, we run the nested sampling simulation.
energies, liveset, _ = nested_sampling(ls, ns_params, 1_000, mc, save)

# ### Metropolis Monte Carlo

# Similarly to nested sampling, we need to define some MC parameters.
# We will use the same Lennard-Jones potential and initial configurations.

# Let's use a walker from the previous example as the initial configuration.
at = ls.walkers[1]

# Define the temperature grid, the number of equilibration steps, the number of sampling steps, and the step size.
temperatures = collect(1000.0:-50:0)
num_equilibration_steps = 100_000
num_sampling_steps = 100_000
step_size = 1.0

# Pass the parameters to the `MetropolisMCParameters` object.
mc_params = MetropolisMCParameters(
    temperatures,
    equilibrium_steps=num_equilibration_steps,
    sampling_steps=num_sampling_steps,
    step_size=step_size,
    step_size_up=1.0,
    accept_range=(0.5,0.5)
)

# Run the Monte Carlo simulation.
mc_energies, mc_ls, mc_cvs, acceptance_rates = monte_carlo_sampling(at, lj, mc_params)


# ### Wang-Landau Sampling

# It's also very easy to set up a Wang-Landau simulation.
# We will use the same Lennard-Jones potential and initial configuration.
at = ls.walkers[1]

# Define the step size and the energy range.
step_size = 1.0
energy_min = -1.26
energy_max = 0.0
num_energy_bins = 1000
# Pass the parameters to the `WangLandauParameters` object.
wl_params = WangLandauParameters(num_steps=10_000,
                            energy_min=energy_min,  
                            energy_max=energy_max, 
                            num_energy_bins=num_energy_bins, 
                            step_size=step_size,
                            max_iter=10000,
                            f_min=1.00001)
# Viola! We can now run the Wang-Landau simulation.
energies_wl, configs, wl_params, S, H = wang_landau(at, lj, wl_params)


# ## Lattice system

# We can also perform simulations on lattice systems.

# Let's set up a 2D square lattice with a supercell of dimensions 4x4x1.
# Make the first four sites occupied.
initial_lattice = SLattice{SquareLattice}(components=[[1,2,3,4]])

# Define the Hamiltonian with the adsorption energy and the nearest-neighbor and next-nearest-neighbor energies.
adsorption_energy = -0.04
nn_energy = -0.01
nnn_energy = -0.0025
h = GenericLatticeHamiltonian(adsorption_energy, [nn_energy, nnn_energy], u"eV")

# ### Wang-Landau Sampling
# This time, let's do a Wang-Landau simulation first, as it's the most straightforward to set up.
energy_min = 20.5 * nn_energy + nn_energy / 8
energy_max = 16 * nn_energy - nn_energy / 8
num_energy_bins = 100

# Pass the parameters to the `WangLandauParameters` object.
wl_params = WangLandauParameters(
    energy_min=energy_min, 
    energy_max=energy_max,
    random_seed=Int(round(time() * 1000)),)

# Run the Wang-Landau simulation.
energies_wl, configs, wl_params, S, H = wang_landau(initial_lattice, h, wl_params)


# ### Metroplis Monte Carlo

# Let's make a Metropolis Monte Carlo simulation on the same lattice system.
mc_lattice = deepcopy(initial_lattice)

# Define the temperature grid, the number of equilibration steps, and the number of sampling steps.
temperatures = collect(200.0:-10:10)  # 200:-1:1
num_equilibration_steps = 25_000
num_sampling_steps = 25_000

# Pass the parameters to the `MetropolisMCParameters` object.
mc_params = MetropolisMCParameters(
    temperatures,
    equilibrium_steps=num_equilibration_steps,
    sampling_steps=num_sampling_steps,
    random_seed=Int(round(time()))
)

# Run the Monte Carlo simulation.
mc_energies, mc_configs, mc_cvs, acceptance_rates = monte_carlo_sampling(mc_lattice, h, mc_params)

# ### Nested Sampling

# Finally, let's set up a nested sampling simulation. It's a bit more complicated than the previous two, 
# and nested sampling is not the optimal method for lattice systems with many degenerate states.

# Let's use the Distributions.jl package to generate the random initial configurations.
using Distributions

# Make 1000 copies of the initial lattice.
walkers = [deepcopy(initial_lattice) for i in 1:1000]

# Now using the `sample` function from the Distributions.jl package, we can randomly select 4 sites to be occupied.
for walker in walkers
    walker.components[1] = [false for i in 1:length(walker.components[1])]
    for i in sample(1:length(walker.components[1]), 4, replace=false)
        walker.components[1][i] = true
    end
end

# Warp the walkers into `LatticeWalker` objects.
walkers = [LatticeWalker(walker) for walker in walkers]

# Again, we use the same Hamiltonian as before.
# We then create a `LatticeGasWalkers` object that contains the walkers and the Hamiltonian.
ls = LatticeGasWalkers(walkers, h)

# The rejection routine for nested sampling is the best choice for lattice systems.
mc = MCRejectionSampling()

# We then set up a `LatticeNestedSamplingParameters` object with the desired parameters.
# We will use the same parameters as before, but we will set the number of Monte Carlo steps to 100.
# We will also set the allowed fail count to 1_000_000, which is a bit high to make sure we explore the configuration space.
# We will also set the random seed to a random number.
ns_params = LatticeNestedSamplingParameters(mc_steps=100,allowed_fail_count=1_000_000,random_seed=Int(floor(1234*rand())))

# Redefine the output file name.
save = SaveEveryN(df_filename="output_ns_lattice2d.csv")

# And finally, we can run the nested sampling simulation
ns_energies, ls, _ = nested_sampling(ls, ns_params, 10_000, mc, save) # src