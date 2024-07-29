# This is a simple example of how to use FreeBird
using FreeBird

np = 8

# First, we need to read the walkers from a file
@elapsed ats = read_walkers("slab.extxyz",pbc="TTF")
@elapsed ats = [AtomWalker{2}(at.configuration;list_num_par=[576,np],frozen=Bool[1,0]) for at in ats]
# Note that the periodic boundary conditions are set to True, True, False

# Then, we need to create a potential object
lj11 = LJParameters(epsilon=0.002844, sigma=3.4, cutoff=3.0, shift=true)
lj22 = LJParameters(epsilon=0.010400, sigma=3.4, cutoff=3.0, shift=true)
lj12 = LJParameters(epsilon=0.005438, sigma=3.4, cutoff=3.0, shift=true)

ljs = CompositeLJParameters(2,[lj11,lj12,lj22])
# Now, we can create the liveset by combining the walkers and the potential
ls = LJAtomWalkers(ats, ljs)
# The last argument is the number of frozen particles

# We can now create a Monte Carlo object
mc = MCRandomWalkClone()

ns_params = NestedSamplingParameters(1600, 0.01, 0.01, 1e-5, 0.1, 0, 200)
# We can also create a save object
save = SaveEveryN(n_snap=10_000)

#slab_height, img_plane, Cs1, Cs2
# And finally, we can run the simulation
@info "starting nested sampling loop"
energies, ls, _ = nested_sampling_loop!(ls, ns_params, 5_000_000, mc, save; slab_height = 3.35u"Å", dispersion_params = [3.35u"Å", 1.6785u"Å", 17.74575u"eV*Å^6", 10.516u"eV*Å^6"])

# The energies are saved in the file output-energies.csv
# We can also save the final walkers to a file
@elapsed write_walkers("slab-$(np)-final-ls.extxyz", ls.walkers)
