using FreeBird
using ExtXYZ

ats = Atoms.(read_frames("scripts/slab.extxyz"))
lj = LennardJonesParameters(cutoff=4.0,shift=true)
liveset = LennardJonesAtomsWalkers(ats, lj)
ns_params = NestedSamplingParameters(100, 0.01, 12)

assign_lj_energies!(liveset)

for i in 1:10000
    emax, _ = nested_sampling_step!(liveset, ns_params)
    println("emax: ", emax)
end