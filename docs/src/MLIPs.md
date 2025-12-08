# Using Machine Learning Interatomic Potentials (MLIPs) in FreeBird

FreeBird.jl supports the use of machine learning interatomic potentials (MLIPs) by interfacing with the ASE (Atomic Simulation Environment) Python calculators. One can use your local Python environment or set up a Conda environment within Julia using `CondaPkg.jl`.

## UMA

For UMA, we recommend setting up a local Conda environment. For example:

```bash
# This is a bash command, run it in your terminal before starting Julia
conda create -prefix ./uma_env python=3.11
conda activate ./uma_env
pip install fairchem-core
```

Now, start Julia with the Conda environment activated, and use `PreferenceTools.jl` to set the `PYTHON` environment variable to point to the Python executable in your Conda environment:

```julia
using PreferenceTools
# press the "]" key to enter Pkg mode
pkg> preference add CondaPkg env=/full/path/to/uma_env
```
Replace `/full/path/to/uma_env` with the actual path to your Conda environment. You only need to do this once.

Now, we can use UMA MLIPs in FreeBird.jl. Here is an example of setting up an `MLIPAtomWalkers` live set:

```julia
using FreeBird
mlp = PyMLPotential(model_type="uma", model_version="uma-s-1p1", task_name="omol", use_gpu=false)
walkers = AtomWalker.(generate_initial_configs(10, 562.5, 6; particle_type=:Si))
ls = MLIPAtomWalkers(walkers, mlp)
```

## MACE

Similarly, for MACE, you can set up a Conda environment and direct CondaPkg.jl to use it. 

Alternatively, you can install the required packages within your Julia environment using CondaPkg.jl:

```julia
using CondaPkg
# press the "]" key to enter Pkg mode
pkg> conda pip_add mace-torch
```

To enable CUDA acceleration with cuEquivariance library, you need to install additional packages:

```julia
pkg> conda pip_add cuequivariance cuequivariance-torch cuequivariance-ops-torch-cu12
```

See [MACE documentation](https://mace.readthedocs.io/en/latest/) for more details on installation and usage.

Now, we can use MACE MLIPs in FreeBird.jl. Here is an example of setting up an `MLIPAtomWalkers` live set:

```julia
using FreeBird
# Create MACE potential using the "small" pre-trained model, you can supply your own model path as well
mlp = PyMLPotential(model_type="mace", model_version="small", float_type="float32", use_gpu=false)
walkers = AtomWalker.(generate_initial_configs(10, 562.5, 6; particle_type=:Si))
ls = MLIPAtomWalkers(walkers, mlp)
```

## ORB
Again, you can point CondaPkg.jl to your local Conda environment where you have installed ORB. 
Alternatively, you can install the required packages within your Julia environment using CondaPkg.jl:

```julia
using CondaPkg
# press the "]" key to enter Pkg mode
pkg> conda pip_add orb-models
```

Now, we can use ORB MLIPs in FreeBird.jl. Here is an example of setting up an `MLIPAtomWalkers` live set:

```julia
using FreeBird
mlp = PyMLPotential(model_type="orb", model_version="n", float_type="float32-high", use_gpu=false)
```

## Add other MLIPs

In principle, you can use any ASE-compatible MLIP by creating a new `PyMLPotential` dispatch in the `FreeBird.AbstractPotentials` module.
Specially, in `src/AbstractPotentials/ase_calculators.jl`, you can add a new function similar to the existing ones for UMA, MACE, and ORB.

Using MACE as an example, you can create a new potential as follows:

```julia
function PyMLPotential(;model_version::String, 
                        float_type::String="",
                        use_gpu::Bool=false)
    # import relevant Python module
    mace = pyimport("mace.calculators")
    # create MACE ASE calculator using an appropriate function call
    mace_calc = mace.mace_mp(model=model_version, default_dtype=float_type, enable_cueq=use_gpu)
    # wrap and return as PyMLPotential
    return PyMLPotential(ASEcalculator(mace_calc))
end
```
FreeBird.jl will then be able to use this new MLIP for evaluating energies, which is defined in the `FreeBird.EnergyEval` module.
Such as:

```julia
function interacting_energy(system::AbstractSystem, calc::PyMLPotential)
    return AtomsCalculators.potential_energy(system, calc.calc)
end
```

## Using MLIPs for nested sampling

Here, we provide an example of using ORB MLIP for nested sampling with `MLIPAtomWalkers` live set.

```julia
using FreeBird

# Set up ORB potential
mlp = PyMLPotential(model_type="orb", model_version="n", float_type="float32-high", use_gpu=false)

# Generate initial configurations and set up walkers, let's say 10 walkers of 6 Si atoms in a box
walkers = AtomWalker.(generate_initial_configs(10, 562.5, 6; particle_type=:Si))

# Set up MLIPAtomWalkers live set, energies will be evaluated using the ORB potential
ls = MLIPAtomWalkers(walkers, mlp)

# Set up nested sampling parameters, using 100 MC steps per iteration and small step sizes
ns_params = NestedSamplingParameters(mc_steps=100, step_size=0.01, step_size_up=0.05, random_seed=1234*rand(Int))

# Set up output saving parameters
save = SaveEveryN(n_traj=10, n_snap=10_000, n_info=1)

# Set up Monte Carlo routines
mc = MCRandomWalkClone() # works with MLIPAtomWalkers on a GPU or CPU; one can use `MCDistributed()` for multi-processing on CPU (no GPU support yet)

# Run nested sampling for 10_000 iterations
energies, liveset, _ = nested_sampling(ls, ns_params, 10_000, mc, save)
```

Fundamentally, using other MLIPs follows the same procedure as above. Be aware of the computational costs with MLIPS, one typically needs to use
GPU for fast energy evaluations, or massively parallel CPU computations to distribute the workload.