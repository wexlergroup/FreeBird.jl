"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV

export read_output
export gamma_factors, partition_function, internal_energy, cv

"""
    read_output(filename::String)

Reads the output file and returns a DataFrame.
"""
function read_output(filename::String)
    return DataFrame(CSV.File(filename))
end

"""
    gamma_factors(iters::Vector{Int}, n_walkers::Int)

Calculates the \$\\Gamma\$ factors for the given number of iterations and walkers.
The \$\\Gamma\$ factors account for the fractions of phase-space volume sampled during
each nested sampling iteration, defined as:
```math
\\Gamma_i = \\frac{1}{N+1} \\left(\\frac{N}{N+1}\\right)^i
```
where \$N\$ is the number of walkers and \$i\$ is the iteration number.

# Arguments
- `iters::Vector{Int}`: The iteration numbers.
- `n_walkers::Int`: The number of walkers.

# Returns
- A vector of \$\\Gamma\$ factors.
"""
function gamma_factors(iters::Vector{Int}, n_walkers::Int)
    gi = (1/(n_walkers+1))*(n_walkers/(n_walkers+1)).^iters
    return gi
end

"""
    partition_function(beta::Float64, gi::Vector{Float64}, ei::Vector{Float64})

Calculates the partition function for the given \$\\beta\$, \$\\Gamma\$ factors, and energies.
The partition function is defined as:
```math
Z(\\beta) = \\sum_i \\Gamma_i \\exp(-E_i \\beta)
```
where \$\\Gamma_i\$ is the \$i\$-th \$\\Gamma\$ factor, \$E_i\$ is the \$i\$-th energy, and \$\\beta\$ is the inverse temperature.

# Arguments
- `beta::Float64`: The inverse temperature.
- `gi::Vector{Float64}`: The \$\\Gamma\$ factors.
- `ei::Vector{Float64}`: The energies.

# Returns
- The partition function.
"""
function partition_function(beta::Float64, 
                            gi::Vector{Float64}, 
                            ei::Vector{Float64})
    z = sum(gi.*exp.(-ei.*beta))
    return z
end

"""
    internal_energy(beta::Float64, gi::Vector{Float64}, ei::Vector{Float64})

Calculates the internal energy from the partition function for the given \$\\beta\$, \$\\Gamma\$ factors, and energies.
The internal energy is defined as:
```math
U(\\beta) = \\frac{\\sum_i \\Gamma_i E_i \\exp(-E_i \\beta)}{\\sum_i \\Gamma_i \\exp(-E_i \\beta)}
```
where \$\\Gamma_i\$ is the \$i\$-th \$\\Gamma\$ factor, \$E_i\$ is the \$i\$-th energy, and \$\\beta\$ is the inverse temperature.

# Arguments
- `beta::Float64`: The inverse temperature.
- `gi::Vector{Float64}`: The \$\\Gamma\$ factors.
- `ei::Vector{Float64}`: The energies.

# Returns
- The internal energy.
"""
function internal_energy(beta::Float64, 
                         gi::Vector{Float64}, 
                         ei::Vector{Float64})
    u = sum(gi.*ei.*exp.(-ei.*beta))/sum(gi.*exp.(-ei.*beta))
    return u
end

"""
    cv(beta::Float64, gi::Vector{Float64}, ei::Vector{Float64}, dof::Int)

Calculates the constant-volume heat capacity at constant volume for the given \$\\beta\$, \$\\Gamma\$ factors, energies, and degrees of freedom.
The heat capacity is defined as:
```math
C_V(\\beta) = \\frac{dof \\cdot k_B}{2} + k_B \\beta^2 \\left(\\frac{\\sum_i \\Gamma_i E_i^2 \\exp(-E_i \\beta)}{Z(\\beta)} - U(\\beta)^2\\right)
```
where \$dof\$ is the degrees of freedom, \$k_B\$ is the Boltzmann constant (in units of eV/K), \$\\beta\$ is the inverse temperature, 
\$\\Gamma_i\$ is the \$i\$-th \$\\Gamma\$ factor, \$E_i\$ is the \$i\$-th energy, \$Z(\\beta)\$ is the partition function, and \$U(\\beta)\$ is the internal energy.

# Arguments
- `beta::Float64`: The inverse temperature.
- `gi::Vector{Float64}`: The \$\\Gamma\$ factors.
- `ei::Vector{Float64}`: The energies.
- `dof::Int`: The degrees of freedom.

# Returns
- The constant-volume heat capacity.
"""
function cv(beta::Float64, 
            gi::Vector{Float64}, 
            ei::Vector{Float64},
            dof::Int)
    z = partition_function(beta, gi, ei)
    u = internal_energy(beta, gi, ei)
    kb = 8.617333262e-5 # eV/K
    cv = dof*kb/2.0 + kb*beta^2 * (sum(gi.*ei.^2 .*exp.(-ei.*beta))/z - u^2)
    return cv
end


end # module AnalysisTools