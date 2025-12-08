"""
    struct PyMLPotential <: SingleComponentPotential{ManyBody}

Light-weight wrapper around an `ase` calculator (accessed via PythonCall) so it
can be passed around in Julia code.

The wrapped calculator *must* be able to return energies in *electron-volts* via
`atoms.get_potential_energy()`.
"""
struct PyMLPotential <: SingleComponentPotential{ManyBody}
    calc::ASEcalculator
end

function PyMLPotential(;model_type::String, 
                        model_version::String, 
                        float_type::String="",
                        task_name::String="",
                        inference_settings::Union{String,Dict}="default",
                        use_gpu::Bool=false)
    os = pyimport("os")
    os.environ["KMP_DUPLICATE_LIB_OK"]=raw"True"

    # Mace
    if model_type == "mace"
        mace = pyimport("mace.calculators")
        mace_calc = mace.mace_mp(model=model_version, default_dtype=float_type, enable_cueq=use_gpu)
        return PyMLPotential(ASEcalculator(mace_calc))
    end

    # Orb
    if model_type == "orb"
        if use_gpu
            device = "cuda"
        else
            device = "cpu"
        end

        pretrained = pyimport("orb_models.forcefield.pretrained")
        ORBCalculator = pyimport("orb_models.forcefield.calculator").ORBCalculator
        orbff = pretrained.orb_v3_direct_20_omat(
            device=device,
            precision=float_type
        )
        orb_calc = ORBCalculator(orbff, device=device)
        return PyMLPotential(ASEcalculator(orb_calc))
    end

    # UMA
    if lowercase(model_type) == "uma" 
        if task_name == ""
            error("Currently, using UMA MLIP model need specify the task_name.")
        end

        if use_gpu
            device = "cuda"
        else
            device = "cpu"
        end

        uma_pretrained = pyimport("fairchem.core.calculate.pretrained_mlip")
        uma_ase_calc = pyimport("fairchem.core.calculate.ase_calculator").FAIRChemCalculator
        
        predictor = uma_pretrained.get_predict_unit(model_name=model_version, 
                                                    inference_settings=inference_settings, 
                                                    device=device)
        uma = uma_ase_calc(predictor, task_name=task_name)
        
        return PyMLPotential(ASEcalculator(uma))
    end
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

