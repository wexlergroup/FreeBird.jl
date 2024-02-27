# This is a simple example of how to use FreeBird
using FreeBird

# First, we need to read the walkers from a file
ats = read_walkers("scripts/slab-08-4x4.extxyz",pbc="TTF")
# Note that the periodic boundary conditions are set to True, True, False

# Then, we need to create a potential object
lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0, shift=0.5)

# Now, we can create the liveset by combining the walkers and the potential
ls = LJAtomWalkers(ats[1:10], lj, 80)
# The last argument is the number of frozen particles

# We can now create a Monte Carlo object
mc = MCDemonWalk(max_add_steps=10000, e_demon_tolerance=0.01u"eV")
mix = MixedMCRoutine(main_routine=MCRandomWalk(),back_up_routine=mc)

# We can also create a save object
save = SaveEveryN("output-energies.csv",100)

# And finally, we can run the simulation
energies, ls, _ = nested_sampling_loop!(ls, 10000, mix, save)

# The energies are saved in the file output-energies.csv
# We can also save the final walkers to a file
write_walkers("slab-08-final-ls.extxyz", ls.walkers)
