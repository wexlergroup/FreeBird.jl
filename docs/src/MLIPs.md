# Using Machine Learning Interatomic Potentials (MLIPs) in FreeBird

## UMA

FreeBird.jl supports the use of machine learning interatomic potentials (MLIPs) by interfacing with the ASE (Atomic Simulation Environment) Python calculators. For UMA, we recommend setting up a local Conda environment. For example:

```bash
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