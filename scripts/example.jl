using FreeBird
using ExtXYZ

ats = Atoms.(read_frames("scripts/slab.extxyz"))
lj = LennardJonesParameters(cutoff=4.0,shift=true)
liveset = LennardJonesAtomsWalkers(ats, lj)
ns_params = NestedSamplingParameters(100, 0.01, 12)

energies, _ = nested_sampling_loop!(liveset, ns_params, 10000)