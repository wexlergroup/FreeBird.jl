# Using Machine Learning Interatomic Potentials (MLIPs) in FreeBird

FreeBird.jl supports the use of machine learning interatomic potentials (MLIPs) by interfacing with the ASE (Atomic Simulation Environment) Python calculators. One can use your local Python environment or set up a Conda environment within Julia using `CondaPkg.jl`.

## FAIRChem UMA model

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
Replace `/full/path/to/uma_env` with the actual path to your Conda environment. You only need to do this once. Restart Julia after setting this preference. To verify that the correct Python environment is being used, you can run:

```julia
using PreferenceTools
# press the "]" key to enter Pkg mode
pkg> preference status
```

You may need to resolve the environment again using:

```julia
using CondaPkg
# press the "]" key to enter Pkg mode
pkg> conda resolve
```

Now, we can use UMA MLIPs in FreeBird.jl. Here is an example of setting up an `MLIPAtomWalkers` live set:

```julia
using FreeBird
mlp = uma_model("omat", "uma-s-1p1", device="cuda")
walkers = AtomWalker.(generate_initial_configs(10, 562.5, 6; particle_type=:Si))
ls = MLIPAtomWalkers(walkers, mlp)
```

The `omat` task name and `uma-s-1p1` model name are required arguments for loading the UMA model. One can pass additional keyword arguments to customize the model loading, for example, `device="cuda"`. See [Fairchem UMA documentation](https://fair-chem.github.io/) for more details on other settings and usage.

## MACE models

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
mlp = mace_model(model_version="small", detype="float32", enable_cueq=true)
```

## ORB models
Again, you can point CondaPkg.jl to your local Conda environment where you have installed ORB. 
Alternatively, you can install the required packages within your Julia environment using CondaPkg.jl:

```julia
using CondaPkg
# press the "]" key to enter Pkg mode
pkg> conda pip_add orb-models
```

Now, we can use ORB MLIPs in FreeBird.jl. Here is an example of setting up an `MLIPAtomWalkers` live set:

```julia
mlp = orb_model(precision="float32-high", device="cuda")
```
See [ORB documentation](https://github.com/orbital-materials/orb-models) for more details on usage.

## CHGNet models

Similarly, you can point CondaPkg.jl to your local Conda environment where you have installed CHGNet.
Alternatively, you can install the required packages within your Julia environment using CondaPkg.jl:
```julia
using CondaPkg
# press the "]" key to enter Pkg mode
pkg> conda pip_add chgnet
```

To use it, simply create a CHGNet model as follows:

```julia
mlp = chgnet_model(model_name="r2scan", use_device="cuda")
```
Here, we use the pre-trained `r2scan` model and enable CUDA. See [CHGNet documentation](https://chgnet.lbl.gov/) for more details on usage.

## UPET models

Again, you can point CondaPkg.jl to the local Conda environment where you have installed UPET.
Alternatively, you can install the required packages within your Julia environment using CondaPkg.jl:
```julia
using CondaPkg
# press the "]" key to enter Pkg mode
pkg> conda pip_add upet
```

We can now load in the UPET model:

```julia
mlp = upet_model(model="pet-mad-s", version="1.0.2", device="cuda")
```
See [UPET documentation](https://github.com/lab-cosmo/upet) for more details on usage.

## Add other MLIPs

In principle, you can use any ASE-compatible MLIP by creating a new `PyMLPotential` dispatch in the `FreeBird.AbstractPotentials` module.
Specially, in `src/AbstractPotentials/ase_calculators.jl`, you can add a new function similar to the existing ones for UMA, MACE, and ORB.

Using MACE as an example, you can create a new potential as follows:

```julia
using FreeBird, PythonCall, ASEconvert
# import relevant Python module
mace = pyimport("mace.calculators")
# create MACE ASE calculator using an appropriate function call
mace_calc = mace.mace_mp(model_version="small", enable_cueq=true)
# wrap and return as PyMLPotential
mlp = PyMLPotential(ASEcalculator(mace_calc))
```
FreeBird.jl will then be able to use this new MLIP for evaluating energies, which is defined in the `FreeBird.EnergyEval` module.
Such as:

```julia
function interacting_energy(system::AbstractSystem, calc::PyMLPotential)
    return AtomsCalculators.potential_energy(system, calc.calc)
end
```

## Using MLIPs for nested sampling

Here, we provide an example of using MACE MLIP for nested sampling with `MLIPAtomWalkers` live set.

```julia
using FreeBird

# Set up MACE potential
mlp = mace_model(model_version="small", detype="float32", enable_cueq=true)

# Generate initial configurations and set up walkers, let's say 10 walkers of 6 Si atoms in a box
walkers = AtomWalker.(generate_initial_configs(10, 562.5, 6; particle_type=:Si))

# Set up MLIPAtomWalkers live set, energies will be evaluated using the ORB potential
ls = MLIPAtomWalkers(walkers, mlp)

# Set up nested sampling parameters, using 100 MC steps per iteration and small step sizes
ns_params = NestedSamplingParameters(mc_steps=100, step_size=0.01, step_size_up=0.05, random_seed=1234*rand(Int))

# Set up output saving parameters
save = SaveEveryN(n_traj=10, n_snap=10_000, n_info=1)

# Set up Monte Carlo routines
mc = MCRandomWalkClone() # works with MLIPAtomWalkers on a GPU or CPU; one can use `MCDistributed()` for multi-processing on CPU (no multi-GPU support yet)

# Run nested sampling for 10_000 iterations
energies, liveset, _ = nested_sampling(ls, ns_params, 10_000, mc, save)
```

Fundamentally, using other MLIPs follows the same procedure as above. Be aware of the computational costs with MLIPS, one typically needs to use GPU for fast energy evaluations, or massively parallel CPU computations to distribute the workload.

## Lattice sampling on an `AtomicLattice`

A `PyMLPotential` can evaluate configurations of an [`AtomicLattice`](@ref).
FreeBird turns each site-occupancy pattern into the corresponding ASE
structure before asking the calculator for its energy.

This example performs a short fixed-particle-number Monte Carlo calculation
for one oxygen atom on a 2×2 Pd(100) surface:

```julia
using FreeBird

mlp = mace_model(
    model="small", default_dtype="float64", enable_cueq=false)

lattice = AtomicLattice{1,SquareLattice}(
    lattice_atom="Pd",
    surface=:fcc100,
    supercell_dimensions=(2, 2, 1),
    lattice_constant=3.947,
    periodicity=(true, true, false),
    adsorbate_atoms=["O"],
    components=[1],
    num_nearest_neighbors=0,
    type_of_sites=["hollow"],
)

params = MetropolisMCParameters(
    [300.0],
    equilibrium_steps=10,
    sampling_steps=20,
    random_seed=42,
)

energies, configurations, heat_capacities, acceptance_rates =
    monte_carlo_sampling(MCNewSample(), lattice, mlp, params)
```

`MCNewSample` moves the adsorbate between lattice sites while keeping its
particle count fixed. Fixed-composition exact enumeration, fixed-N nested
sampling, and Wang–Landau sampling can use the same lattice and potential.
Multi-species lattices are also supported when the underlying calculator can
evaluate every listed element.

The third entry of `supercell_dimensions` is the number of substrate layers,
not a vacuum thickness. Converge the layer count for the property of interest,
and validate or relax `adsorbate_height` for the selected substrate, site,
adsorbate, and potential. The constructor applies one height to every species
and does not perform structural relaxation.

The example uses double precision because Monte Carlo decisions depend on
energy differences between configurations. Single precision can be faster, but
should be treated as an explicit performance tradeoff and checked against a
double-precision run.

MLIP lattice moves evaluate the complete atomic structure for every changed
proposal, so calculations can be much slower than a lattice Hamiltonian.
During each temperature run, the sampler stores an independent ASE frame at
every step, so peak memory grows with both slab size and step count.
Start with a small surface and short run when checking a new model. For lattice
sampling, use `MLattice` with a `ClassicalHamiltonian` when no
explicit atomic structure is needed; a `PyMLPotential` lattice calculation
requires an `AtomicLattice`.
