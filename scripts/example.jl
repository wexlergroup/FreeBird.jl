# This is a simple example of how to use FreeBird
using FreeBird

# First, we need to read the walkers from a file
ats = read_walkers("scripts/slab.extxyz",pbc="TTF")
# Note that the periodic boundary conditions are set to True, True, False

# Then, we need to create a potential object
lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0, shift=true)

# Now, we can create the liveset by combining the walkers and the potential
ls = LJAtomWalkers(ats, lj, 12)
# The last argument is the number of frozen particles

# We can now create a Monte Carlo object
mc = MCRandomWalkClone()

# Setting the parameters for the nested sampling
ns_params = NestedSamplingParameters(1600, 0.1, 0.01, 1e-5, 1.0, 0, 200)

# We can also create a save object
save = SaveEveryN(n=1000)

# And finally, we can run the simulation (takes about 1 minute)
energies, ls, _ = nested_sampling_loop!(ls, ns_params, 10_000, mc, save)

# The energies are saved in the file output-energies.csv
# We can also save the final walkers to a file
write_walkers("output-final-ls.extxyz", ls.walkers)
