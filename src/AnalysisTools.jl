"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV

export read_output
export ωᵢ, partition_function, internal_energy, cv

"""
    read_output(filename::String)

Reads the output file and returns a DataFrame.
"""
function read_output(filename::String)
    return DataFrame(CSV.File(filename))
end

"""
    ωᵢ(iters::Vector{Int}, n_walkers::Int)

Calculates the \$\\omega\$ factors for the given number of iterations and walkers.
The \$\\omega\$ factors account for the fractions of parameter-space volume sampled during
each nested sampling iteration, defined as:
```math
\\omega_i = \\frac{1}{N+1} \\left(\\frac{N}{N+1}\\right)^i
```
where \$N\$ is the number of walkers and \$i\$ is the iteration number.

# Arguments
- `iters::Vector{Int}`: The iteration numbers.
- `n_walkers::Int`: The number of walkers.

# Returns
- A vector of \$\\omega\$ factors.
"""
function ωᵢ(iters::Vector{Int}, n_walkers::Int)
    ω_i = (1/(n_walkers+1))*(n_walkers/(n_walkers+1)).^iters
    return ω_i
end

"""
    partition_function(β::Float64, ωi::Vector{Float64}, Ei::Vector{Float64})

Calculates the partition function for the given \$\\beta\$, \$\\omega\$ factors, and energies.
The partition function is defined as:
```math
Z(\\beta) = \\sum_i \\omega_i \\exp(-E_i \\beta)
```
where \$\\omega_i\$ is the \$i\$-th \$\\omega\$ factor, \$E_i\$ is the \$i\$-th energy, and \$\\beta\$ is the inverse temperature.

# Arguments
- `β::Float64`: The inverse temperature.
- `ωi::Vector{Float64}`: The \$\\omega\$ factors.
- `Ei::Vector{Float64}`: The energies.

# Returns
- The partition function.
"""
function partition_function(β::Float64, 
                            ωi::Vector{Float64}, 
                            Ei::Vector{Float64})
    z = sum(ωi.*exp.(-Ei.*β))
    return z
end

"""
    internal_energy(β::Float64, ωi::Vector{Float64}, ei::Vector{Float64})

Calculates the internal energy from the partition function for the given \$\\beta\$, \$\\omega\$ factors, and energies.
The internal energy is defined as:
```math
U(\\beta) = \\frac{\\sum_i \\omega_i E_i \\exp(-E_i \\beta)}{\\sum_i \\omega_i \\exp(-E_i \\beta)}
```
where \$\\omega_i\$ is the \$i\$-th \$\\omega\$ factor, \$E_i\$ is the \$i\$-th energy, and \$\\beta\$ is the inverse temperature.

# Arguments
- `β::Float64`: The inverse temperature.
- `ωi::Vector{Float64}`: The \$\\omega\$ factors.
- `Ei::Vector{Float64}`: The energies in eV.

# Returns
- The internal energy.
"""
function internal_energy(β::Float64, 
                         ωi::Vector{Float64}, 
                         Ei::Vector{Float64})
    u = sum(ωi.*Ei.*exp.(-Ei.*β))/sum(ωi.*exp.(-Ei.*β))
    return u
end

"""
    cv(β::Float64, omega_i::Vector{Float64}, Ei::Vector{Float64}, dof::Int)

Calculates the constant-volume heat capacity at constant volume for the given \$\\beta\$, \$\\omega\$ factors, energies, and degrees of freedom.
The heat capacity is defined as:
```math
C_V(\\beta) = \\frac{dof \\cdot k_B}{2} + k_B \\beta^2 \\left(\\frac{\\sum_i \\omega_i E_i^2 \\exp(-E_i \\beta)}{Z(\\beta)} - U(\\beta)^2\\right)
```
where \$dof\$ is the degrees of freedom, \$k_B\$ is the Boltzmann constant (in units of eV/K), \$\\beta\$ is the inverse temperature, 
\$\\omega_i\$ is the \$i\$-th \$\\omega\$ factor, \$E_i\$ is the \$i\$-th energy, \$Z(\\beta)\$ is the partition function, and \$U(\\beta)\$ is the internal energy.

# Arguments
- `β::Float64`: The inverse temperature.
- `ωi::Vector{Float64}`: The \$\\omega\$ factors.
- `Ei::Vector{Float64}`: The energies in eV.
- `dof::Int`: The degrees of freedom.

# Returns
- The constant-volume heat capacity.
"""
function cv(β::Float64, 
            ωi::Vector{Float64}, 
            Ei::Vector{Float64},
            dof::Int)
    z = partition_function(β, ωi, Ei)
    u = internal_energy(β, ωi, Ei)
    kb = 8.617333262e-5 # eV/K
    cv = dof*kb/2.0 + kb*β^2 * (sum(ωi.*Ei.^2 .*exp.(-Ei.*β))/z - u^2)
    return cv
end

"""
    cv(df::DataFrame, βs::Vector{Float64}, dof::Int, n_walkers::Int)

Calculates the constant-volume heat capacity at constant volume for the given DataFrame, inverse temperatures, degrees of freedom, and number of walkers.
The heat capacity is defined as:
```math
C_V(\\beta) = \\frac{dof \\cdot k_B}{2} + k_B \\beta^2 \\left(\\frac{\\sum_i \\omega_i E_i^2 \\exp(-E_i \\beta)}{Z(\\beta)} - U(\\beta)^2\\right)
```
where \$dof\$ is the degrees of freedom, \$k_B\$ is the Boltzmann constant (in units of eV/K), \$\\beta\$ is the inverse temperature,
\$\\omega_i\$ is the \$i\$-th \$\\omega\$ factor, \$E_i\$ is the \$i\$-th energy, \$Z(\\beta)\$ is the partition function, and \$U(\\beta)\$ is the internal energy.

# Arguments
- `df::DataFrame`: The DataFrame containing the output data.
- `βs::Vector{Float64}`: The inverse temperatures.
- `dof::Int`: The degrees of freedom.
- `n_walkers::Int`: The number of walkers.

# Returns
- A vector of constant-volume heat capacities.
"""
function cv(df::DataFrame, βs::Vector{Float64}, dof::Int, n_walkers::Int)
    ωi = ωᵢ(df.iter, n_walkers)
    Ei = df.emax .- minimum(df.emax)
    cvs = [cv(b, ωi, Ei, dof) for b in βs]
    return cvs
end


end # module AnalysisTools