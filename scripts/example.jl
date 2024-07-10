# This is a simple example of how to use FreeBird
using FreeBird

# First, we need to read the walkers from a file
# Note that the periodic boundary conditions are updated to True, True, False
ats = read_walkers("scripts/slab-08-4x4.extxyz",pbc="TTF")

# Update the system to have 80 frozen particles and 8 free particles
ats = [AtomWalker{2}(at.configuration;list_num_par=[80,8],frozen=Bool[1,0]) for at in ats]

# Then, we need to create a potential object
lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0, shift=true)

# Now, we can create the liveset by combining the walkers and the potential
ls = LJAtomWalkers(ats, lj)

# We can now create a Monte Carlo object
mc = MCRandomWalkClone()

ns_params = NestedSamplingParameters(1600, 0.1, 0.01, 1e-5, 0.1, 0, 200)
# We can also create a save object, using default value for most parameters except n_snap
save = SaveEveryN(n_snap=1000)

# And finally, we can run the simulation
energies, ls, _ = nested_sampling_loop!(ls, ns_params, 10_000, mc, save)

# Save the final walkers to a file
write_walkers("slab-08-final-ls.extxyz", ls.walkers)
