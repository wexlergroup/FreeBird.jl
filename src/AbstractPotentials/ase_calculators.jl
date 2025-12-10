"""
    struct PyMLPotential <: SingleComponentPotential{ManyBody}

Light-weight wrapper around an `ase` calculator (accessed via PythonCall) so it
can be passed around in Julia code.

The wrapped calculator *must* be able to return energies in *electron-volts* via
`atoms.get_potential_energy()`.

# Available Models
- MACE: see [`mace_model`](@ref) function)
- ORB: see [`orb_model`](@ref) function)
- UMA: see [`uma_model`](@ref) function)
- CHGNet: see [`chgnet_model`](@ref) function)
"""
struct PyMLPotential <: SingleComponentPotential{ManyBody}
    calc::ASEcalculator
end

"""
    mace_model(arg...; kwargs...) <: PyMLPotential

Loading MACE foundation model Python module through ASE interface.
See https://github.com/ACEsuit/mace-foundations for the available arguments and usage.

# Example
mlp = mace_model(model_version="small", detype="float32", enable_cueq=true)
"""
function mace_model(arg...; kwargs...)
    mace = pyimport("mace.calculators")
    mace_calc = mace.mace_mp(arg...; kwargs...)
    return PyMLPotential(ASEcalculator(mace_calc))
end

"""
    orb_model(arg...; kwargs...) <: PyMLPotential

Loading ORB potential Python module through ASE interface.
See https://github.com/orbital-materials/orb-models for the available arguments and usage.

# Example
mlp = orb_model(precision="float32-high", device="cuda")
"""
function orb_model(arg...; kwargs...)
    pretrained = pyimport("orb_models.forcefield.pretrained")
    ORBCalculator = pyimport("orb_models.forcefield.calculator").ORBCalculator
    orbff = pretrained.orb_v3_direct_20_omat(arg...; kwargs...)
    orb_calc = ORBCalculator(orbff)
    return PyMLPotential(ASEcalculator(orb_calc))
end

"""
    uma_model(task_name::String, model_name::String; kwargs...) <: PyMLPotential

Loading Fairchem UMA model Python module through ASE interface.
See https://fair-chem.github.io/ for the available arguments and usage.
`task_name` and `model_name` are required.

# Example
mlp = uma_model("omat", "uma-s-1p1", device="cuda")
"""
function uma_model(task_name::String, model_name::String; kwargs...)
    uma_pretrained = pyimport("fairchem.core.calculate.pretrained_mlip")
    uma_ase_calc = pyimport("fairchem.core.calculate.ase_calculator").FAIRChemCalculator

    predictor = uma_pretrained.get_predict_unit(model_name=model_name; kwargs...)
    uma = uma_ase_calc(predictor, task_name=task_name)
    return PyMLPotential(ASEcalculator(uma))
end

"""
    chgnet_model(arg...; kwargs...) <: PyMLPotential

Loading CHGNet model Python module through ASE interface.
See https://chgnet.lbl.gov/ for the available arguments and usage.

# Example
mlp = chgnet_model(model_name="r2scan", use_device="cuda")
"""
function chgnet_model(arg...; kwargs...)
    model = pyimport("chgnet.model.model")
    pot = model.CHGNet.load(arg...; kwargs...)
    calc = pyimport("chgnet.model.dynamics").CHGNetCalculator(potential=pot, property="energy")
    return PyMLPotential(ASEcalculator(calc))
end

