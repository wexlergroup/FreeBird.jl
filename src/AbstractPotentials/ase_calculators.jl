"""
    struct ASECalculator

Light-weight wrapper around an `ase` calculator (accessed via PythonCall) so it
can be passed around in Julia code just like the built-in `LJParameters`.

The wrapped calculator *must* be able to return energies in *electron-volts* via
`atoms.get_potential_energy()`.
"""
struct PyCalculator <: SingleComponentPotential{ManyBody}
    calc::ASEcalculator
end

"""
    ASELennardJones(;epsilon=1.0, sigma=1.0, cutoff=Inf)

Convenience constructor that creates an `ase.calculators.lj.LennardJones`
calculator with the supplied parameters and wraps it in an `ASECalculator`.
All energy/length units follow the ASE convention (eV / Å).
"""
function ASELennardJones(; epsilon = 1.0, sigma = 1.0, cutoff = Inf)
    lj_mod = pyimport("ase.calculators.lj")
    # ASE uses the keyword `rc` for the cut-off radius
    rc = isfinite(cutoff) ? cutoff * sigma : cutoff
    py_calc = lj_mod.LennardJones(; epsilon = epsilon, sigma = sigma, rc = rc)
    return ASECalculator(py_calc)
end


function MACEPotential(; mace_model::String, default_dtype::String="float64")
    os = pyimport("os")
    os.environ["KMP_DUPLICATE_LIB_OK"]="True"
    mace = pyimport("mace.calculators")
    mace_calc = mace.mace_mp(model=mace_model, default_dtype=default_dtype)
    return ASEcalculator(mace_calc)
end

function OrbPotential(; device::String="cpu", precision::String="float64")
    os = pyimport("os")
    os.environ["KMP_DUPLICATE_LIB_OK"]="True"
    pretrained = pyimport("orb_models.forcefield.pretrained")
    ORBCalculator = pyimport("orb_models.forcefield.calculator.ORBCalculator")
    orbff = pretrained.orb_v3_conservative_inf_omat(
        device=device,
        precision=precision
    )
    orb_calc = ORBCalculator(orbff, device=device)
    return ASEcalculator(orb_calc)
end

# Dummy frozen_energy for compatibility with LJAtomWalkers constructor when no
# atoms are actually frozen (returns zero).
function frozen_energy(system::AtomsBase.AbstractSystem,
                       calc::PyCalculator,
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool})
    # For now we simply return 0 because ASE calculators do not currently
    # support separating frozen–frozen interactions.
    return 0.0u"eV"
end

