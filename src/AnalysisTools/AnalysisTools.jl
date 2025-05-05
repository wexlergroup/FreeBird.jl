"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV, Arrow

export read_output
export ωᵢ, partition_function, internal_energy, cv

"""
    read_output(filename::String)

Reads the output file and returns a DataFrame.
"""
function read_output(filename::String)
    if splitext(filename)[end] == ".csv"
        data = CSV.File(filename)
    elseif splitext(filename)[end] == ".arrow"
        data = Arrow.Table(filename)
    else
        error("Unsupported file format. Please provide a .csv or .arrow file.")
    end
    return DataFrame(data)
end

"""
    ωᵢ(iters::Vector{Int}, n_walkers::Int; n_cull::Int=1, ω0::Float64=1.0)

Calculates the \$\\omega\$ factors for the given number of iterations and walkers.
The \$\\omega\$ factors account for the fractions of phase-space volume sampled during
each nested sampling iteration, defined as:
```math
\\omega_i = \\frac{C}{K+C} \\left(\\frac{K}{K+C}\\right)^i
```
where \$K\$ is the number of walkers, \$C\$ is the number of culled walkers, 
and \$i\$ is the iteration number.

# Arguments
- `iters::Vector{Int}`: The iteration numbers.
- `n_walkers::Int`: The number of walkers.
- `n_cull::Int`: The number of culled walkers. Default is 1.
- `ω0::Float64`: The initial \$\\omega\$ factor. Default is 1.0.

# Returns
- A vector of \$\\omega\$ factors.
"""
function ωᵢ(iters::Vector{Int}, n_walkers::Int; n_cull::Int=1, ω0::Float64=1.0)
    ωi = ω0 * (n_cull/(n_walkers+n_cull)) * (n_walkers/(n_walkers+n_cull)).^iters
    return ωi
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

Calculates the constant-volume heat capacity for the given \$\\beta\$, \$\\omega\$ factors, energies, and degrees of freedom.
The heat capacity is defined as:
```math
C_V(\\beta) = \\frac{\\mathrm{dof} \\cdot k_B}{2} + k_B \\beta^2 \\left(\\frac{\\sum_i \\omega_i E_i^2 \\exp(-E_i \\beta)}{Z(\\beta)} - U(\\beta)^2\\right)
```
where \$\\mathrm{dof}\$ is the degrees of freedom, \$k_B\$ is the Boltzmann constant (in units of eV/K), \$\\beta\$ is the inverse temperature, 
\$\\omega_i\$ is the \$i\$-th \$\\omega\$ factor, \$E_i\$ is the \$i\$-th energy, \$Z(\\beta)\$ is the partition function, and \$U(\\beta)\$ is the internal energy.

# Arguments
- `β::Float64`: The inverse temperature.
- `ωi::Vector{Float64}`: The \$\\omega\$ factors.
- `Ei::Vector{Float64}`: The energies in eV.
- `dof::Int`: The degrees of freedom, equals to the number of dimensions times the number of particles.

# Returns
- The constant-volume heat capacity.
"""
function cv(β::Float64,
            ωi::Vector{Float64}, 
            Ei::Vector{Float64},
            dof::Int64;
            kb::Float64=8.617333262e-5)
    expo = ωi.*exp.(-Ei.*β)
    ei_expo = Ei.*expo
    ei2_expo = Ei.*ei_expo
    z = sum(expo)
    u = sum(ei_expo)/z
    cv = dof*kb/2.0 + kb*β^2 * (sum(ei2_expo)/z - u^2)
    return cv
end

"""
    cv(df::DataFrame, βs::Vector{Float64}, dof::Int, n_walkers::Int)

(Nested Sampling) Calculates the constant-volume heat capacity at constant volume for the given DataFrame, inverse temperatures, degrees of freedom, and number of walkers.
The heat capacity is defined as:
```math
C_V(\\beta) = \\frac{\\mathrm{dof} \\cdot k_B}{2} + k_B \\beta^2 \\left(\\frac{\\sum_i \\omega_i E_i^2 \\exp(-E_i \\beta)}{Z(\\beta)} - U(\\beta)^2\\right)
```
where \$\\mathrm{dof}\$ is the degrees of freedom, \$k_B\$ is the Boltzmann constant (in units of eV/K), \$\\beta\$ is the inverse temperature,
\$\\omega_i\$ is the \$i\$-th \$\\omega\$ factor, \$E_i\$ is the \$i\$-th energy, \$Z(\\beta)\$ is the partition function, and \$U(\\beta)\$ is the internal energy.

# Arguments
- `df::DataFrame`: The DataFrame containing the output data.
- `βs::Vector{Float64}`: The inverse temperatures.
- `dof::Int`: The degrees of freedom, equals to the number of dimensions times the number of particles. For a lattice, it is zero.
- `n_walkers::Int`: The number of walkers.
- `n_cull::Int`: The number of culled walkers. Default is 1.
- `ω0::Float64`: The initial \$\\omega\$ factor. Default is 1.0.

# Returns
- A vector of constant-volume heat capacities.
"""
function cv(df::DataFrame, 
            βs::Vector{Float64}, 
            dof::Int, n_walkers::Int; 
            n_cull::Int=1, 
            ω0::Float64=1.0, 
            kb::Float64=8.617333262e-5)
    ωi = ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
    Ei = df.emax .- minimum(df.emax)
    cvs = Vector{Float64}(undef, length(βs))
    Threads.@threads for (i, b) in collect(enumerate(βs))
        cvs[i] = cv(b, ωi, Ei, dof; kb=kb)
    end
    return cvs
end

"""
    cv(Ts::Vector{Float64}, dof::Int, energy_bins::Vector{Float64}, entropy::Vector{Float64})

(Wang-Landau Sampling) Calculates the constant-volume heat capacity at constant volume for the given temperatures, degrees of 
freedom, energy bins, and entropy. The kinetic energy is treated classically, and is added to the heat capacity as \$dof \\cdot k_B/2\$.

# Arguments
- `Ts::Vector{Float64}`: The temperatures in Kelvin.
- `dof::Int`: The degrees of freedom, equals to the number of dimensions times the number of particles. For a lattice, it is zero.
- `energy_bins::Vector{Float64}`: The energy bins in eV.
- `entropy::Vector{Float64}`: The entropy.

# Returns
- A vector of constant-volume heat capacities.
"""
function cv(Ts::Vector{Float64}, dof::Int, energy_bins::Vector{Float64}, entropy::Vector{Float64})
    kb = 8.617333262e-5 # eV/K
    β = 1 ./(kb.*Ts)
    E_rel = energy_bins .- minimum(energy_bins)
    S_shifted = entropy .- minimum(entropy[entropy .> 0])
    g = exp.(S_shifted)
    Z = zeros(length(Ts))
    E_avg = zeros(length(Ts))
    E2_avg = zeros(length(Ts))
    Cv = zeros(length(Ts))
    for (i, temp) in enumerate(Ts)
        Z[i] = sum(exp.(-E_rel ./ (kb * temp)) .* g)
        E_avg[i] = sum(E_rel .* exp.(-E_rel ./ (kb * temp)) .* g) / Z[i]
        E2_avg[i] = sum(E_rel.^2 .* exp.(-E_rel ./ (kb * temp)) .* g) / Z[i]
        Cv[i] = (E2_avg[i] .- E_avg[i].^2) ./ (kb * temp.^2) .+ dof*kb/2
    end
    return Cv
end


end # module AnalysisTools