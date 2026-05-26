"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV, Arrow
using Statistics: mean

export read_output
export ωᵢ, partition_function, internal_energy, cv
export gc_thermodynamic_stats
export ideal_gas_reference_thermodynamic_stats
export effective_sample_size_ideal_gas_reference
export combined_ideal_gas_reference_thermodynamic_stats

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
function ωᵢ(iters::AbstractVector{Int}, n_walkers::Int; n_cull::Int=1, ω0::Float64=1.0)
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
- `βs::Vector{Float64}`: Inverse temperatures.
- `dof::Int`: The degrees of freedom, equal to the number of dimensions times the number of particles. For a lattice, it is zero.
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


"""
    gc_thermodynamic_stats(β::Float64, ωi::Vector{Float64},
                           grand_energies::Vector{Float64},
                           energies::Vector{Float64},
                           numbers::Vector{Int},
                           μ::Float64;
                           kb::Float64=8.617333262e-5)

Compute grand-canonical thermodynamic averages from nested sampling output.

The log-sum-exp trick is used for numerical stability. The grand-canonical
heat capacity at constant μ is:

    C_{V,μ} = k_B β² [Var(E) − μ Cov(E, N)]

# Arguments
- `β::Float64`: Inverse temperature 1/(k_B T).
- `ωi::Vector{Float64}`: Phase-space volume weights from NS.
- `grand_energies::Vector{Float64}`: Ω_i = E_i − μ N_i values.
- `energies::Vector{Float64}`: E_i values.
- `numbers::Vector{Int}`: N_i values.
- `μ::Float64`: Chemical potential.
- `kb::Float64`: Boltzmann constant (default: eV/K).

# Returns
- `(⟨E⟩, C_{V,μ}, ⟨N⟩)`: Mean energy, GC heat capacity, mean particle number.
"""
function gc_thermodynamic_stats(β::Float64,
                                 ωi::Vector{Float64},
                                 grand_energies::Vector{Float64},
                                 energies::Vector{Float64},
                                 numbers::Vector{Int},
                                 μ::Float64;
                                 kb::Float64=8.617333262e-5)
    n = length(ωi)
    if n != length(grand_energies) || n != length(energies) || n != length(numbers)
        throw(DimensionMismatch("All input vectors must have the same length"))
    end
    if n == 0
        return NaN, NaN, NaN
    end

    # Log-sum-exp for numerical stability
    log_terms = [log(ωi[i]) - β * grand_energies[i] for i in 1:n]
    max_log = maximum(log_terms)

    z = 0.0
    u = 0.0   # ⟨E⟩
    u2 = 0.0  # ⟨E²⟩
    n_sum = 0.0  # ⟨N⟩
    en_sum = 0.0 # ⟨EN⟩

    for i in 1:n
        w = exp(log_terms[i] - max_log)
        z += w
        u += w * energies[i]
        u2 += w * energies[i]^2
        n_sum += w * numbers[i]
        en_sum += w * energies[i] * numbers[i]
    end

    if z == 0.0
        return NaN, NaN, NaN
    end

    u /= z
    u2 /= z
    n_avg = n_sum / z
    en_avg = en_sum / z

    var_e = u2 - u^2
    cov_en = en_avg - u * n_avg

    cv = kb * β^2 * (var_e - μ * cov_en)

    return u, cv, n_avg
end

"""
    gc_thermodynamic_stats(df::DataFrame, βs::Vector{Float64},
                           n_walkers::Int, μ::Float64;
                           n_cull::Int=1, ω0::Float64=1.0,
                           kb::Float64=8.617333262e-5)

Compute grand-canonical thermodynamic stats from a GC-NS output DataFrame.

The DataFrame must have columns `:iter`, `:omega`, `:energy`, `:num_particles`.
Adds the live-walker contribution to the end of the recorded samples for correct
normalization.

# Arguments
- `df::DataFrame`: GC-NS output with columns `[:iter, :omega, :energy, :num_particles]`.
- `βs::Vector{Float64}`: Inverse temperatures at which to evaluate.
- `n_walkers::Int`: Number of walkers used in the NS run.
- `μ::Float64`: Chemical potential.
- `n_cull::Int=1`: Number of walkers culled per iteration.
- `ω0::Float64=1.0`: Initial phase-space volume.
- `kb::Float64`: Boltzmann constant (default: eV/K).

# Returns
- `(mean_E, Cv, mean_N)`: Vectors of ⟨E⟩, C_{V,μ}, and ⟨N⟩ at each β.
"""
function gc_thermodynamic_stats(df::DataFrame,
                                 βs::Vector{Float64},
                                 n_walkers::Int,
                                 μ::Float64;
                                 n_cull::Int=1,
                                 ω0::Float64=1.0,
                                 kb::Float64=8.617333262e-5)
    ωi = ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
    grand_es = df.omega
    Es = df.energy
    Ns = df.num_particles

    mean_Es = Vector{Float64}(undef, length(βs))
    Cvs = Vector{Float64}(undef, length(βs))
    mean_Ns = Vector{Float64}(undef, length(βs))

    Threads.@threads for (i, b) in collect(enumerate(βs))
        mean_Es[i], Cvs[i], mean_Ns[i] = gc_thermodynamic_stats(
            b, ωi, grand_es, Es, Ns, μ; kb=kb)
    end

    return mean_Es, Cvs, mean_Ns
end


# ----------------------------------------------------------------------
# Internal log-weight helpers for ideal-gas-reference reweighting.
# Shared by ideal_gas_reference_thermodynamic_stats,
# effective_sample_size_ideal_gas_reference, and the multi-run
# combination functions.
# ----------------------------------------------------------------------

"""
    _ideal_gas_reference_log_boltzmann_weights(β, μ, ωi, U, N; z0)

Return the per-sample log Boltzmann weights for reweighting a U-sorted
ideal-gas-reference NS run to target ``(\\mu, T)``:

```math
\\log \\tilde w_j = \\log \\omega_j + N_j (\\beta \\mu - \\log z_0) - \\beta U_j.
```

Internal helper. Validates input lengths and `z0 > 0`.
"""
function _ideal_gas_reference_log_boltzmann_weights(β::Float64, μ::Float64,
                                                     ωi::Vector{Float64},
                                                     U::Vector{Float64},
                                                     N::AbstractVector{<:Integer};
                                                     z0::Float64)
    n = length(ωi)
    if n != length(U) || n != length(N)
        throw(DimensionMismatch("All input vectors must have the same length"))
    end
    if z0 <= 0.0
        throw(ArgumentError("z0 must be positive, got $z0"))
    end

    log_z0 = log(z0)
    log_terms = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        log_terms[i] = log(ωi[i]) + N[i] * (β * μ - log_z0) - β * U[i]
    end
    return log_terms
end

"""
    _ideal_gas_reference_log_importance_ratio(β, μ, U, N; z0, T, T0=nothing, kb=…)

Return the per-sample log importance ratio (without the NS phase-space
factor ``\\omega``) for reweighting from ``(z_0, T_0)`` to target
``(\\mu, T)``:

```math
\\log r_j(\\mu, T; T_0) = N_j (\\log z(\\mu,T) - \\log z_0)
                       - U_j \\left(\\tfrac{1}{k_B T} - \\tfrac{1}{k_B T_0}\\right).
```

When `T0 === nothing`, defaults to `T0 = T` (no temperature reweighting,
only the fugacity ratio). Internal helper.
"""
function _ideal_gas_reference_log_importance_ratio(β::Float64, μ::Float64,
                                                    U::Vector{Float64},
                                                    N::AbstractVector{<:Integer};
                                                    z0::Float64,
                                                    T::Float64,
                                                    T0::Union{Float64,Nothing}=nothing,
                                                    kb::Float64=8.617333262e-5)
    n = length(U)
    if n != length(N)
        throw(DimensionMismatch("U and N must have the same length"))
    end
    if z0 <= 0.0
        throw(ArgumentError("z0 must be positive, got $z0"))
    end
    if T <= 0.0
        throw(ArgumentError("T must be positive, got $T"))
    end
    if T0 !== nothing && T0 <= 0.0
        throw(ArgumentError("T0 must be positive, got $T0"))
    end

    T0_eff = T0 === nothing ? T : T0
    log_z0 = log(z0)
    log_z = β * μ
    inv_kbT = 1.0 / (kb * T)
    inv_kbT0 = 1.0 / (kb * T0_eff)
    delta_inv_kbT = inv_kbT - inv_kbT0

    log_terms = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        log_terms[i] = N[i] * (log_z - log_z0) - U[i] * delta_inv_kbT
    end
    return log_terms
end


"""
    ideal_gas_reference_thermodynamic_stats(β::Float64, μ::Float64,
                                             ωi::Vector{Float64},
                                             U::Vector{Float64},
                                             N::AbstractVector{<:Integer};
                                             z0::Float64,
                                             kb::Float64=8.617333262e-5,
                                             return_variances::Bool=false)

Compute grand-canonical thermodynamic averages from a U-sorted,
ideal-gas-reference nested sampling output, reweighted to a target
``(\\mu, T)``.

The samples were drawn against the prior ``\\pi_0(\\sigma) \\propto z_0^{N(\\sigma)}``.
Reweighting to the target chemical potential ``\\mu`` and inverse temperature
``\\beta = 1/(k_B T)`` uses importance weights ``(z/z_0)^{N_j}`` with
``z = e^{\\beta\\mu}``. The log-Boltzmann weight per sample is

```math
\\log \\tilde w_j = \\log \\omega_j + N_j (\\beta \\mu - \\log z_0) - \\beta U_j.
```

Sums are evaluated via the log-sum-exp trick for numerical stability.
The returned ``\\ln Z_{\\mathrm{relative}}`` is the log of the
NS-estimated partition sum *relative to the prior normalization*; the
absolute grand partition function on a lattice with ``M`` sites is
``\\ln \\Xi(\\mu, T) = M \\ln(1 + z_0) + \\ln Z_{\\mathrm{relative}}``,
generalizing the ``M \\ln 2`` offset of the Ω-sorted GCNS (which
corresponds to ``z_0 = 1``).

The grand-canonical heat capacity at constant ``\\mu`` is

```math
C_{V,\\mu} = k_B \\beta^2 \\big(\\mathrm{Var}(U) - \\mu \\,\\mathrm{Cov}(U, N)\\big).
```

# Arguments
- `β::Float64`: Inverse temperature ``1/(k_B T)``.
- `μ::Float64`: Target chemical potential.
- `ωi::Vector{Float64}`: NS phase-space-volume weights per recorded sample.
- `U::Vector{Float64}`: Potential energies per recorded sample.
- `N::AbstractVector{<:Integer}`: Particle counts per recorded sample.
- `z0::Float64`: Reference fugacity used during sampling.
- `kb::Float64`: Boltzmann constant (default eV/K).
- `return_variances::Bool=false`: If `true`, additionally return per-target
  importance-sampling variances for `mean_E`, `mean_N`, and
  `ln_Z_relative` (see "Returns" below). Variances for `Cv_mu` are not
  provided.

# Returns
- If `return_variances=false` (default):
  `(mean_E, mean_N, Cv_mu, ln_Z_relative)::NTuple{4,Float64}`.
- If `return_variances=true`:
  `(mean_E, mean_N, Cv_mu, ln_Z_relative,
    var_mean_E, var_mean_N, var_ln_Z_relative)::NTuple{7,Float64}`.

The variances follow standard importance-sampling delta-method approximations:
- `var_mean_E ≈ Var_w(U) / N_eff`
- `var_mean_N ≈ Var_w(N) / N_eff`
- `var_ln_Z_relative ≈ 1/N_eff − 1/N_total`

where ``N_{\\mathrm{eff}} = (\\sum_j \\tilde w_j)^2 / \\sum_j \\tilde w_j^2``
is the importance-sampling effective sample size at ``(\\mu, T)``. These
estimators are accurate when the importance weights are near-uniform; in
extreme reweighting regimes (very small ``N_{\\mathrm{eff}}``) they
under-estimate the true uncertainty.
"""
function ideal_gas_reference_thermodynamic_stats(β::Float64, μ::Float64,
                                                  ωi::Vector{Float64},
                                                  U::Vector{Float64},
                                                  N::AbstractVector{<:Integer};
                                                  z0::Float64,
                                                  kb::Float64=8.617333262e-5,
                                                  return_variances::Bool=false)
    n = length(ωi)
    if n != length(U) || n != length(N)
        throw(DimensionMismatch("All input vectors must have the same length"))
    end
    if z0 <= 0.0
        throw(ArgumentError("z0 must be positive, got $z0"))
    end
    if n == 0
        return return_variances ?
               (NaN, NaN, NaN, NaN, NaN, NaN, NaN) :
               (NaN, NaN, NaN, NaN)
    end

    log_terms = _ideal_gas_reference_log_boltzmann_weights(β, μ, ωi, U, N; z0=z0)
    max_log = maximum(log_terms)

    z = 0.0
    z2_shifted = 0.0  # Σ exp(2 (log_terms - max_log)) for N_eff
    e_sum = 0.0
    e2_sum = 0.0
    n_sum = 0.0
    n2_sum = 0.0
    en_sum = 0.0

    @inbounds for i in 1:n
        shifted = log_terms[i] - max_log
        w = exp(shifted)
        z += w
        z2_shifted += exp(2 * shifted)
        e_sum += w * U[i]
        e2_sum += w * U[i]^2
        n_sum += w * N[i]
        n2_sum += w * N[i]^2
        en_sum += w * U[i] * N[i]
    end

    if z == 0.0 || !isfinite(z)
        return return_variances ?
               (NaN, NaN, NaN, NaN, NaN, NaN, NaN) :
               (NaN, NaN, NaN, NaN)
    end

    e_avg = e_sum / z
    e2_avg = e2_sum / z
    n_avg = n_sum / z
    n2_avg = n2_sum / z
    en_avg = en_sum / z

    var_e = e2_avg - e_avg^2
    var_n = n2_avg - n_avg^2
    cov_en = en_avg - e_avg * n_avg

    Cv_mu = kb * β^2 * (var_e - μ * cov_en)
    ln_Z_relative = max_log + log(z)

    if !return_variances
        return e_avg, n_avg, Cv_mu, ln_Z_relative
    end

    # Importance-sampling effective sample size at this (μ, T).
    # N_eff = (Σ w)^2 / Σ w^2; the max_log shift cancels in the ratio.
    inv_N_eff = z2_shifted / (z * z)

    # Clamp slightly-negative roundoff variances to zero before scaling.
    var_mean_E = max(var_e, 0.0) * inv_N_eff
    var_mean_N = max(var_n, 0.0) * inv_N_eff
    # Delta-method: Var(ln Σ w) ≈ 1/N_eff − 1/N_total.
    var_ln_Z = max(0.0, inv_N_eff - 1.0 / n)

    return e_avg, n_avg, Cv_mu, ln_Z_relative, var_mean_E, var_mean_N, var_ln_Z
end

"""
    ideal_gas_reference_thermodynamic_stats(df::DataFrame,
                                             μs::AbstractVector{<:Real},
                                             Ts::AbstractVector{<:Real},
                                             n_walkers::Int, z0::Float64;
                                             n_cull::Int=1, ω0::Float64=1.0,
                                             kb::Float64=8.617333262e-5,
                                             return_variances::Bool=false)

Sweep [`ideal_gas_reference_thermodynamic_stats`](@ref) over a 2D grid
of ``(\\mu, T)`` from a U-sorted GCNS output DataFrame.

The DataFrame must have columns `:iter`, `:emax`, and `:num_particles`.

# Returns
- If `return_variances=false` (default): a NamedTuple with fields
  `mean_E`, `mean_N`, `Cv_mu`, `ln_Z_relative` — each a
  `length(μs) × length(Ts)` matrix.
- If `return_variances=true`: the above plus `var_mean_E`, `var_mean_N`,
  `var_ln_Z_relative` matrices (delta-method importance-sampling
  variances; see scalar method for details and limitations).
"""
function ideal_gas_reference_thermodynamic_stats(df::DataFrame,
                                                  μs::AbstractVector{<:Real},
                                                  Ts::AbstractVector{<:Real},
                                                  n_walkers::Int, z0::Float64;
                                                  n_cull::Int=1, ω0::Float64=1.0,
                                                  kb::Float64=8.617333262e-5,
                                                  return_variances::Bool=false)
    ωi = ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
    Us = Vector{Float64}(df.emax)
    Ns = Vector{Int}(df.num_particles)

    nμ = length(μs)
    nT = length(Ts)
    mean_E = Matrix{Float64}(undef, nμ, nT)
    mean_N = Matrix{Float64}(undef, nμ, nT)
    Cv_mu = Matrix{Float64}(undef, nμ, nT)
    ln_Z_rel = Matrix{Float64}(undef, nμ, nT)

    if return_variances
        var_E = Matrix{Float64}(undef, nμ, nT)
        var_N = Matrix{Float64}(undef, nμ, nT)
        var_lnZ = Matrix{Float64}(undef, nμ, nT)
    end

    for j in 1:nT
        β = 1.0 / (kb * Ts[j])
        for i in 1:nμ
            μ = Float64(μs[i])
            if return_variances
                (mean_E[i,j], mean_N[i,j], Cv_mu[i,j], ln_Z_rel[i,j],
                 var_E[i,j], var_N[i,j], var_lnZ[i,j]) =
                    ideal_gas_reference_thermodynamic_stats(
                        β, μ, ωi, Us, Ns; z0=z0, kb=kb,
                        return_variances=true)
            else
                (mean_E[i,j], mean_N[i,j], Cv_mu[i,j], ln_Z_rel[i,j]) =
                    ideal_gas_reference_thermodynamic_stats(
                        β, μ, ωi, Us, Ns; z0=z0, kb=kb)
            end
        end
    end

    if return_variances
        return (mean_E=mean_E, mean_N=mean_N, Cv_mu=Cv_mu, ln_Z_relative=ln_Z_rel,
                var_mean_E=var_E, var_mean_N=var_N, var_ln_Z_relative=var_lnZ)
    else
        return (mean_E=mean_E, mean_N=mean_N, Cv_mu=Cv_mu, ln_Z_relative=ln_Z_rel)
    end
end

"""
    effective_sample_size_ideal_gas_reference(β::Float64, μ::Float64,
                                               ωi::Vector{Float64},
                                               U::Vector{Float64},
                                               N::AbstractVector{<:Integer};
                                               z0::Float64,
                                               T::Float64,
                                               T0::Union{Float64,Nothing}=nothing,
                                               kb::Float64=8.617333262e-5)

Compute the importance-sampling effective sample size for reweighting a
U-sorted ideal-gas-reference NS run from its sampling distribution
``(z_0, T_0)`` to a target ``(\\mu, T)``.

Per Eq. 3 of the companion proposal,

```math
r_j(\\mu, T; T_0) = (z(\\mu,T)/z_0)^{N_j} \\exp\\!\\left[-\\tfrac{U_j}{k_B}\\!\\left(\\tfrac{1}{T}-\\tfrac{1}{T_0}\\right)\\right],
\\quad
N_{\\mathrm{eff}}(\\mu, T) = \\frac{(\\sum_j \\omega_j r_j)^2}{\\sum_j \\omega_j^2 r_j^2}.
```

When `T0 === nothing`, defaults to `T0 = T` (the importance weight then
contains only the fugacity ratio ``(z/z_0)^N``, with no temperature
reweighting). Computation uses log-sum-exp on the log-weights for
numerical stability.

# Arguments
- `β::Float64`: Target inverse temperature ``1/(k_B T)``.
- `μ::Float64`: Target chemical potential.
- `ωi::Vector{Float64}`: NS phase-space-volume weights per recorded sample.
- `U::Vector{Float64}`: Potential energies per recorded sample.
- `N::AbstractVector{<:Integer}`: Particle counts per recorded sample.
- `z0::Float64`: Reference fugacity used during sampling.
- `T::Float64`: Target temperature (K), used only via `1/T - 1/T0`. Should
  satisfy ``\\beta = 1/(k_B T)`` consistently with the first argument.
- `T0::Union{Float64,Nothing}`: Sampling reference temperature. `nothing`
  ⇒ defaults to `T`.
- `kb::Float64`: Boltzmann constant (default eV/K).

# Returns
- `N_eff::Float64`. Returns `NaN` if the sample is empty or if denominators vanish.
"""
function effective_sample_size_ideal_gas_reference(β::Float64, μ::Float64,
                                                    ωi::Vector{Float64},
                                                    U::Vector{Float64},
                                                    N::AbstractVector{<:Integer};
                                                    z0::Float64,
                                                    T::Float64,
                                                    T0::Union{Float64,Nothing}=nothing,
                                                    kb::Float64=8.617333262e-5)
    n = length(ωi)
    if n != length(U) || n != length(N)
        throw(DimensionMismatch("All input vectors must have the same length"))
    end
    if n == 0
        return NaN
    end

    # Helper validates z0, T, T0 and returns log r_j (without ω).
    log_r = _ideal_gas_reference_log_importance_ratio(
        β, μ, U, N; z0=z0, T=T, T0=T0, kb=kb)

    # Σ_j ω_j r_j and Σ_j ω_j^2 r_j^2 via log-sum-exp.
    log_num_terms = Vector{Float64}(undef, n)
    log_den_terms = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        log_num_terms[i] = log(ωi[i]) + log_r[i]
        log_den_terms[i] = 2 * log(ωi[i]) + 2 * log_r[i]
    end

    max_num = maximum(log_num_terms)
    max_den = maximum(log_den_terms)

    s_num = 0.0
    s_den = 0.0
    @inbounds for i in 1:n
        s_num += exp(log_num_terms[i] - max_num)
        s_den += exp(log_den_terms[i] - max_den)
    end

    if s_num == 0.0 || s_den == 0.0 || !isfinite(s_num) || !isfinite(s_den)
        return NaN
    end

    log_neff = 2 * (max_num + log(s_num)) - (max_den + log(s_den))
    return exp(log_neff)
end

"""
    effective_sample_size_ideal_gas_reference(df::DataFrame,
                                               μs::AbstractVector{<:Real},
                                               Ts::AbstractVector{<:Real},
                                               n_walkers::Int, z0::Float64;
                                               T0::Union{Float64,Nothing}=nothing,
                                               n_cull::Int=1, ω0::Float64=1.0,
                                               kb::Float64=8.617333262e-5)

Sweep [`effective_sample_size_ideal_gas_reference`](@ref) over a 2D grid
of ``(\\mu, T)`` from a U-sorted GCNS output DataFrame.

# Returns
- `N_eff::Matrix{Float64}` of size `length(μs) × length(Ts)`.
"""
function effective_sample_size_ideal_gas_reference(df::DataFrame,
                                                    μs::AbstractVector{<:Real},
                                                    Ts::AbstractVector{<:Real},
                                                    n_walkers::Int, z0::Float64;
                                                    T0::Union{Float64,Nothing}=nothing,
                                                    n_cull::Int=1, ω0::Float64=1.0,
                                                    kb::Float64=8.617333262e-5)
    ωi = ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
    Us = Vector{Float64}(df.emax)
    Ns = Vector{Int}(df.num_particles)

    nμ = length(μs)
    nT = length(Ts)
    Neff = Matrix{Float64}(undef, nμ, nT)

    for j in 1:nT
        T = Float64(Ts[j])
        β = 1.0 / (kb * T)
        for i in 1:nμ
            μ = Float64(μs[i])
            Neff[i,j] = effective_sample_size_ideal_gas_reference(
                β, μ, ωi, Us, Ns; z0=z0, T=T, T0=T0, kb=kb)
        end
    end

    return Neff
end


# ----------------------------------------------------------------------
# Multi-run combination: pool K U-sorted ideal-gas-reference NS runs
# at distinct reference fugacities z₀^(k) into a single (μ, T)-resolved
# estimate.
#
# Two methods:
#   :inverse_variance — combine per-run estimates by 1/σ² weighting
#   :mixture_importance — pooled mixture-importance sampling (Veach-Guibas
#                        balance heuristic adapted for NS, with full
#                        absolute-Ξ correction)
#
# Both consume the per-run output of
# `ideal_gas_reference_grand_canonical_nested_sampling`.
# ----------------------------------------------------------------------

"""
    _logsumexp(xs::AbstractVector{Float64}) -> Float64

Numerically stable ``\\log \\sum_i \\exp(x_i)``.

Returns `-Inf` for an empty input. If the maximum element is `+Inf` or
`-Inf`, the result is the maximum (no normalization possible). `NaN`
elements propagate.

Internal helper.
"""
function _logsumexp(xs::AbstractVector{Float64})
    isempty(xs) && return -Inf
    m = maximum(xs)
    isfinite(m) || return m
    s = 0.0
    @inbounds for x in xs
        s += exp(x - m)
    end
    return m + log(s)
end

"""
    _log_mixture_prior_density(N_j::Integer,
                                log_alpha_over_total::Vector{Float64},
                                log_z0s::Vector{Float64}) -> Float64

Compute ``\\log f_{\\mathrm{mix}}^{\\mathrm{unn}}(\\sigma_j)
= \\log \\sum_k (\\alpha_k/N_{\\mathrm{total}}) \\, z_0^{(k)\\,N_j}``,
i.e. the (unnormalized) mixture-prior log-density for a sample with
particle count `N_j`, evaluated against `K` runs with per-run sample
counts `α_k` and reference fugacities `z₀^{(k)}`.

The two input vectors are precomputed once per combination call:
- `log_alpha_over_total[k] = log(α_k / N_total)`
- `log_z0s[k]              = log(z₀^{(k)})`

Computed via log-sum-exp for numerical stability. Allocation-free
(O(K) work per call). Internal helper.
"""
function _log_mixture_prior_density(N_j::Integer,
                                     log_alpha_over_total::Vector{Float64},
                                     log_z0s::Vector{Float64})
    K = length(log_alpha_over_total)
    @assert K == length(log_z0s) "log_alpha_over_total and log_z0s must have the same length"
    @assert K ≥ 1 "at least one run is required"

    # Find max term first (manual LSE to avoid allocating a Vector).
    max_term = log_alpha_over_total[1] + N_j * log_z0s[1]
    @inbounds for k in 2:K
        t = log_alpha_over_total[k] + N_j * log_z0s[k]
        if t > max_term
            max_term = t
        end
    end

    s = 0.0
    @inbounds for k in 1:K
        t = log_alpha_over_total[k] + N_j * log_z0s[k]
        s += exp(t - max_term)
    end

    return max_term + log(s)
end


"""
    _inverse_variance_combine_one(x_k, var_k, valid)

Inverse-variance combine a per-run scalar quantity. Returns the tuple
`(combined, var_combined)`. Runs with `valid[k]=false`, non-finite `x_k`,
or non-finite negative `var_k` are excluded.

If at least one surviving run has `var_k = 0` exactly (e.g., the
non-interacting Tier-0 case where `Var(U) = 0`), it dominates with
infinite precision: combined value is the equal-weight mean of the
zero-variance runs, combined variance is `0`.

Returns `(NaN, NaN)` if no run survives. Internal helper.
"""
function _inverse_variance_combine_one(x_k::Vector{Float64},
                                        var_k::Vector{Float64},
                                        valid::AbstractVector{Bool})
    K = length(x_k)
    keep = falses(K)
    @inbounds for k in 1:K
        keep[k] = valid[k] && isfinite(x_k[k]) && isfinite(var_k[k]) && var_k[k] ≥ 0.0
    end
    if !any(keep)
        return NaN, NaN
    end

    # Special case: degenerate σ²=0 runs — they dominate.
    zero_var_mask = keep .& (var_k .== 0.0)
    if any(zero_var_mask)
        idx = findall(zero_var_mask)
        return mean(x_k[idx]), 0.0
    end

    # Standard inverse-variance.
    idx = findall(keep)
    precisions = 1.0 ./ var_k[idx]
    precision_sum = sum(precisions)
    if !isfinite(precision_sum) || precision_sum == 0.0
        return NaN, NaN
    end
    weighted_sum = sum(x_k[idx] .* precisions)
    return weighted_sum / precision_sum, 1.0 / precision_sum
end

"""
    _gather_per_run_stats(dfs, z0s, μs, Ts, n_walkers; n_cull, ω0, kb)

Compute per-run thermodynamic stats (with variances) and per-run N_eff
matrices for each of the K runs. Returns
`(per_run_stats::Vector{NamedTuple}, per_run_neff::Array{Float64,3})`,
where `per_run_neff[i, j, k]` is `N_eff^(k)` at `(μs[i], Ts[j])`.

Internal helper used by `combined_ideal_gas_reference_thermodynamic_stats`.
"""
function _gather_per_run_stats(dfs::Vector{DataFrame},
                                z0s::Vector{Float64},
                                μs::AbstractVector{<:Real},
                                Ts::AbstractVector{<:Real},
                                n_walkers::Int;
                                n_cull::Int=1, ω0::Float64=1.0,
                                kb::Float64=8.617333262e-5)
    K = length(dfs)
    nμ = length(μs)
    nT = length(Ts)
    per_run_stats = Vector{NamedTuple}(undef, K)
    per_run_neff = Array{Float64,3}(undef, nμ, nT, K)
    for k in 1:K
        per_run_stats[k] = ideal_gas_reference_thermodynamic_stats(
            dfs[k], μs, Ts, n_walkers, z0s[k];
            n_cull=n_cull, ω0=ω0, kb=kb, return_variances=true)
        per_run_neff[:, :, k] = effective_sample_size_ideal_gas_reference(
            dfs[k], μs, Ts, n_walkers, z0s[k];
            n_cull=n_cull, ω0=ω0, kb=kb)
    end
    return per_run_stats, per_run_neff
end

"""
    _validate_combined_inputs(dfs, z0s, n_sites, min_N_eff, method)

Validate inputs to `combined_ideal_gas_reference_thermodynamic_stats`.
Throws `ArgumentError` on length mismatch, non-positive `z0`, etc.
Internal helper.
"""
function _validate_combined_inputs(dfs::Vector{DataFrame},
                                    z0s::Vector{Float64},
                                    n_sites::Int,
                                    min_N_eff::Float64,
                                    method::Symbol)
    K = length(dfs)
    if K != length(z0s)
        throw(ArgumentError("length(dfs) ($K) must equal length(z0s) ($(length(z0s)))"))
    end
    if K < 1
        throw(ArgumentError("at least one run is required"))
    end
    for (k, z) in enumerate(z0s)
        if z ≤ 0.0
            throw(ArgumentError("z0s[$k] = $z must be positive"))
        end
    end
    if n_sites < 1
        throw(ArgumentError("n_sites must be at least 1, got $n_sites"))
    end
    if min_N_eff < 0.0
        throw(ArgumentError("min_N_eff must be non-negative, got $min_N_eff"))
    end
    if method ∉ (:inverse_variance, :mixture_importance)
        throw(ArgumentError("method must be :inverse_variance or :mixture_importance, got :$method"))
    end
    return nothing
end


"""
    combined_ideal_gas_reference_thermodynamic_stats(
        dfs::Vector{DataFrame}, z0s::Vector{Float64},
        μs::AbstractVector{<:Real}, Ts::AbstractVector{<:Real},
        n_walkers::Int, n_sites::Int;
        method::Symbol=:mixture_importance,
        min_N_eff::Float64=0.0,
        n_cull::Int=1, ω0::Float64=1.0,
        kb::Float64=8.617333262e-5,
        return_per_run::Bool=false)

Combine `K` U-sorted ideal-gas-reference GCNS runs at distinct reference
fugacities into a single ``(\\mu, T)``-resolved estimate of the grand-
canonical thermodynamic observables.

Each input DataFrame must have columns `[:iter, :emax, :num_particles]`,
matching the schema produced by
[`ideal_gas_reference_grand_canonical_nested_sampling`](@ref). The
fugacities `z0s[k]` must be paired with the corresponding `dfs[k]`.

# Methods

- `:inverse_variance` — Combine per-run estimates with weights ``1/\\sigma_k^2``,
  where ``\\sigma_k^2`` is the per-run importance-sampling variance at
  ``(\\mu, T)``. Conservative and reduces exactly to single-run estimation
  at `K = 1`. The combined `ln_Z_relative` mixes per-run additive
  constants `M·ln(1+z₀^(k))` and is therefore not directly convertible
  to absolute Ξ — use `:mixture_importance` for that.
- `:mixture_importance` — Pooled multiple-importance-sampling (Veach-Guibas
  balance heuristic) against the mixture prior
  ``\\pi_{\\mathrm{mix}}(\\sigma) = \\sum_k (\\alpha_k/N_{\\mathrm{total}}) \\pi_0^{(k)}(\\sigma)``.
  Returns absolute ``\\ln \\Xi`` via the mixture-aware additive constant
  ``\\mathrm{LSE}_k[\\log(\\alpha_k/N_{\\mathrm{total}}) + M \\log(1 + z_0^{(k)})]``.
  *(This method is implemented in step 4 of the campaign; calling now
  raises `ErrorException`.)*

# Per-run N_eff filtering

When `min_N_eff > 0`, runs with `N_eff^(k)(μ, T) < min_N_eff` are
excluded *per (μ, T) cell* (not per whole run). For
`:inverse_variance`, excluded runs drop from the inverse-variance sum.
If every run is excluded at some cell, the output for that cell is
`NaN`. Cell-by-cell filtering is a regularization; document if you use
it.

# Returns

A `NamedTuple`. Always present:
- `mean_E`, `mean_N`, `Cv_mu`, `ln_Z_relative` — `length(μs) × length(Ts)`
  matrices (μ index first, T second).

For `:inverse_variance`:
- `var_mean_E`, `var_mean_N`, `var_ln_Z_relative` — combined variances.
- `N_eff_per_run :: Array{Float64,3}` of size `(nμ, nT, K)` — per-run N_eff.
- `Cv_mu` is the simple mean of valid runs' `Cv_mu` (no variance is
  produced; document as a diagnostic only).

If `return_per_run=true`:
- `per_run :: Vector{NamedTuple}` — per-run single-run estimates
  (output of `ideal_gas_reference_thermodynamic_stats(df_k, ...; return_variances=true)`).

# 1-run reduction

With `K = 1`, `:inverse_variance` reduces exactly to single-run
estimation:
`combined_..._stats([df], [z0], μs, Ts, K_walk, M; method=:inverse_variance)`
matches `ideal_gas_reference_thermodynamic_stats(df, μs, Ts, K_walk, z0; return_variances=true)`
to within floating-point tolerance.
"""
function combined_ideal_gas_reference_thermodynamic_stats(
        dfs::Vector{DataFrame},
        z0s::Vector{Float64},
        μs::AbstractVector{<:Real},
        Ts::AbstractVector{<:Real},
        n_walkers::Int,
        n_sites::Int;
        method::Symbol=:mixture_importance,
        min_N_eff::Float64=0.0,
        n_cull::Int=1,
        ω0::Float64=1.0,
        kb::Float64=8.617333262e-5,
        return_per_run::Bool=false)
    _validate_combined_inputs(dfs, z0s, n_sites, min_N_eff, method)

    per_run_stats, per_run_neff = _gather_per_run_stats(
        dfs, z0s, μs, Ts, n_walkers; n_cull=n_cull, ω0=ω0, kb=kb)

    if method === :inverse_variance
        return _combine_inverse_variance(
            per_run_stats, per_run_neff, μs, Ts;
            min_N_eff=min_N_eff, return_per_run=return_per_run)
    elseif method === :mixture_importance
        return _combine_mixture_importance(
            dfs, z0s, μs, Ts, n_walkers, n_sites;
            min_N_eff=min_N_eff, n_cull=n_cull, ω0=ω0, kb=kb,
            return_per_run=return_per_run,
            per_run_stats=per_run_stats, per_run_neff=per_run_neff)
    end
end

"""
    _combine_inverse_variance(per_run_stats, per_run_neff, μs, Ts;
                              min_N_eff, return_per_run)

Method A combination over a 2D (μ, T) grid. Internal helper.
"""
function _combine_inverse_variance(per_run_stats::Vector{<:NamedTuple},
                                    per_run_neff::Array{Float64,3},
                                    μs::AbstractVector{<:Real},
                                    Ts::AbstractVector{<:Real};
                                    min_N_eff::Float64,
                                    return_per_run::Bool)
    K = length(per_run_stats)
    nμ = length(μs)
    nT = length(Ts)

    mean_E_combined = Matrix{Float64}(undef, nμ, nT)
    mean_N_combined = Matrix{Float64}(undef, nμ, nT)
    Cv_mu_combined = Matrix{Float64}(undef, nμ, nT)
    ln_Z_combined = Matrix{Float64}(undef, nμ, nT)
    var_E_combined = Matrix{Float64}(undef, nμ, nT)
    var_N_combined = Matrix{Float64}(undef, nμ, nT)
    var_lnZ_combined = Matrix{Float64}(undef, nμ, nT)

    e_buf = Vector{Float64}(undef, K)
    n_buf = Vector{Float64}(undef, K)
    cv_buf = Vector{Float64}(undef, K)
    lnz_buf = Vector{Float64}(undef, K)
    var_e_buf = Vector{Float64}(undef, K)
    var_n_buf = Vector{Float64}(undef, K)
    var_lnz_buf = Vector{Float64}(undef, K)
    valid = Vector{Bool}(undef, K)

    for j in 1:nT, i in 1:nμ
        for k in 1:K
            stats = per_run_stats[k]
            e_buf[k] = stats.mean_E[i, j]
            n_buf[k] = stats.mean_N[i, j]
            cv_buf[k] = stats.Cv_mu[i, j]
            lnz_buf[k] = stats.ln_Z_relative[i, j]
            var_e_buf[k] = stats.var_mean_E[i, j]
            var_n_buf[k] = stats.var_mean_N[i, j]
            var_lnz_buf[k] = stats.var_ln_Z_relative[i, j]
            neff_k = per_run_neff[i, j, k]
            valid[k] = isfinite(neff_k) && neff_k ≥ min_N_eff
        end

        mean_E_combined[i, j], var_E_combined[i, j] =
            _inverse_variance_combine_one(e_buf, var_e_buf, valid)
        mean_N_combined[i, j], var_N_combined[i, j] =
            _inverse_variance_combine_one(n_buf, var_n_buf, valid)
        ln_Z_combined[i, j], var_lnZ_combined[i, j] =
            _inverse_variance_combine_one(lnz_buf, var_lnz_buf, valid)

        # Cv_mu: simple mean of valid runs' Cv (no variance available).
        valid_cv_idx = findall(k -> valid[k] && isfinite(cv_buf[k]), 1:K)
        Cv_mu_combined[i, j] = isempty(valid_cv_idx) ? NaN : mean(cv_buf[valid_cv_idx])
    end

    base = (mean_E=mean_E_combined,
            mean_N=mean_N_combined,
            Cv_mu=Cv_mu_combined,
            ln_Z_relative=ln_Z_combined,
            var_mean_E=var_E_combined,
            var_mean_N=var_N_combined,
            var_ln_Z_relative=var_lnZ_combined,
            N_eff_per_run=per_run_neff)
    return return_per_run ? merge(base, (per_run=per_run_stats,)) : base
end


"""
    _combine_mixture_importance(dfs, z0s, μs, Ts, n_walkers, n_sites; ...)

Method B combination over a 2D (μ, T) grid via the Veach-Guibas balance
heuristic against the **normalized** mixture prior

```math
\\pi_{\\mathrm{mix}}^{\\mathrm{norm}}(\\sigma)
   = \\sum_k (\\alpha_k/N_{\\mathrm{total}}) \\pi_0^{(k)}(\\sigma)
   = \\sum_k (\\alpha_k/N_{\\mathrm{total}}) \\, z_0^{(k)\\,N(\\sigma)} / (1+z_0^{(k)})^M.
```

The pooled NS weights ``\\tilde\\omega_j := (\\alpha_{k_j}/N_{\\mathrm{total}}) \\omega_j``
satisfy ``\\sum_j \\tilde\\omega_j \\,g(\\sigma_j) \\approx \\langle g \\rangle_{\\pi_{\\mathrm{mix}}^{\\mathrm{norm}}}``
for any test function ``g``. Setting ``g = f_{\\mathrm{target}} / \\pi_{\\mathrm{mix}}^{\\mathrm{norm}}``
gives an unbiased estimator of ``\\Xi(\\mu, T)`` directly, with
**no extra additive constant**:

```math
\\log \\tilde w_j = \\log(\\alpha_{k_j}/N_{\\mathrm{total}}) + \\log \\omega_j
                  + N_j \\log z(\\mu, T) - \\beta U_j
                  - \\log \\pi_{\\mathrm{mix}}^{\\mathrm{norm}}(\\sigma_j),
\\quad
\\log \\pi_{\\mathrm{mix}}^{\\mathrm{norm}}(\\sigma_j) =
   \\mathrm{LSE}_k\\!\\left[\\log(\\alpha_k/N_{\\mathrm{total}}) - M \\log(1+z_0^{(k)}) + N_j \\log z_0^{(k)}\\right].
```

A naive `f_{\\mathrm{mix}}^{\\mathrm{unn}}`-based importance weight (without
the per-run ``M \\log(1+z_0^{(k)})`` term inside the LSE) is **biased for
K≥2 with distinct z₀**: the ratio ``\\pi_{\\mathrm{mix}}^{\\mathrm{norm}}/f_{\\mathrm{mix}}^{\\mathrm{unn}}``
is not σ-independent, so no single additive constant recovers Ξ.

Reductions per cell:
- ``\\langle A \\rangle = \\sum_j \\tilde w_j A_j / \\sum_j \\tilde w_j``
- ``\\ln \\Xi_{\\mathrm{absolute}} = \\mathrm{LSE}_j(\\log \\tilde w_j)`` (direct, no additive constant)
- ``\\ln \\Xi_{\\mathrm{relative}} = \\ln \\Xi_{\\mathrm{absolute}} - \\log Z_{\\mathrm{mix}}``
  with ``\\log Z_{\\mathrm{mix}} = \\mathrm{LSE}_k\\!\\left[\\log(\\alpha_k/N_{\\mathrm{total}}) + M \\log(1+z_0^{(k)})\\right]``
  (bookkeeping; reduces to ``M \\log(1+z_0)`` at K=1, matching the single-run convention)
- ``N_{\\mathrm{eff,pooled}} = (\\sum_j \\tilde w_j)^2 / \\sum_j \\tilde w_j^2``

When per-cell `min_N_eff` filtering excludes runs, the mixture is recomputed
over the surviving subset only (samples from excluded runs are dropped
from the pool, and ``N_{\\mathrm{total}}`` is recomputed). Otherwise the
"all-runs" mixture quantities are precomputed once and reused.

Internal helper.
"""
function _combine_mixture_importance(dfs::Vector{DataFrame},
                                      z0s::Vector{Float64},
                                      μs::AbstractVector{<:Real},
                                      Ts::AbstractVector{<:Real},
                                      n_walkers::Int,
                                      n_sites::Int;
                                      min_N_eff::Float64,
                                      n_cull::Int,
                                      ω0::Float64,
                                      kb::Float64,
                                      return_per_run::Bool,
                                      per_run_stats::Vector{<:NamedTuple},
                                      per_run_neff::Array{Float64,3})
    K = length(dfs)
    nμ = length(μs)
    nT = length(Ts)
    M = n_sites

    # Per-run sample-level data (avoid recomputation in the cell loop).
    ωi_per_run = [ωᵢ(dfs[k].iter, n_walkers; n_cull=n_cull, ω0=ω0) for k in 1:K]
    Us_per_run = [Vector{Float64}(dfs[k].emax) for k in 1:K]
    Ns_per_run = [Vector{Int}(dfs[k].num_particles) for k in 1:K]
    log_ω_per_run = [log.(ωi_per_run[k]) for k in 1:K]

    # Per-run sample counts and run-level mixture quantities.
    αs = Float64[length(Ns_per_run[k]) for k in 1:K]
    N_total = sum(αs)
    log_z0s = log.(z0s)
    log_alpha_all = log.(αs ./ N_total)
    # Normalized-mixture per-run weight: log(α_k/N_total) − M·log(1+z₀^(k)).
    # Passing this to `_log_mixture_prior_density` makes its LSE compute
    # log π_mix^norm(σ) instead of log f_mix^unn(σ); see derivation in the
    # docstring above.
    log_normed_weights_all = [log_alpha_all[k] - M * log1p(z0s[k]) for k in 1:K]

    # log Z_mix is now bookkeeping only — used to convert ln_Z_absolute to
    # the single-run-style ln_Z_relative. The pooled estimator is unbiased
    # for ln_Z_absolute directly (no additive constant required).
    log_Z_mix_all_terms = [log_alpha_all[k] + M * log1p(z0s[k]) for k in 1:K]
    log_Z_mix_all = _logsumexp(log_Z_mix_all_terms)

    # All-runs per-sample log π_mix^norm(σ_j) — depends only on N_j and
    # the run-level mixture quantities, so precompute once and reuse
    # across cells when no per-cell filter excludes any run.
    log_pi_mix_all = Vector{Vector{Float64}}(undef, K)
    for k in 1:K
        Ns_k = Ns_per_run[k]
        v = Vector{Float64}(undef, length(Ns_k))
        for j in 1:length(Ns_k)
            v[j] = _log_mixture_prior_density(Ns_k[j], log_normed_weights_all, log_z0s)
        end
        log_pi_mix_all[k] = v
    end

    # Output matrices.
    mean_E = Matrix{Float64}(undef, nμ, nT)
    mean_N = Matrix{Float64}(undef, nμ, nT)
    Cv_mu = Matrix{Float64}(undef, nμ, nT)
    ln_Z_rel = Matrix{Float64}(undef, nμ, nT)
    ln_Z_abs = Matrix{Float64}(undef, nμ, nT)
    N_eff_pooled = Matrix{Float64}(undef, nμ, nT)

    for j_T in 1:nT, i_μ in 1:nμ
        T = Float64(Ts[j_T])
        β = 1.0 / (kb * T)
        μ_target = Float64(μs[i_μ])

        # Per-cell valid mask
        valid = falses(K)
        for k in 1:K
            neff_k = per_run_neff[i_μ, j_T, k]
            valid[k] = isfinite(neff_k) && neff_k ≥ min_N_eff
        end

        if !any(valid)
            mean_E[i_μ, j_T]   = NaN
            mean_N[i_μ, j_T]   = NaN
            Cv_mu[i_μ, j_T]    = NaN
            ln_Z_rel[i_μ, j_T] = NaN
            ln_Z_abs[i_μ, j_T] = NaN
            N_eff_pooled[i_μ, j_T] = NaN
            continue
        end

        # Use precomputed mixture if no run filtered; otherwise rebuild
        # for the valid subset.
        if all(valid)
            valid_run_indices = collect(1:K)
            log_alpha_eff = log_alpha_all
            log_Z_mix_eff = log_Z_mix_all
            log_pi_mix_eff = log_pi_mix_all
        else
            valid_run_indices = findall(valid)
            αs_valid = αs[valid_run_indices]
            N_total_valid = sum(αs_valid)
            log_alpha_eff = log.(αs_valid ./ N_total_valid)
            log_z0s_valid = log_z0s[valid_run_indices]
            # See derivation in the docstring: pass the normalized-mixture
            # weight `log(α/N_total) − M·log(1+z₀)` to the helper so its LSE
            # produces log π_mix^norm.
            log_normed_weights_eff = [
                log_alpha_eff[vk] - M * log1p(z0s[valid_run_indices[vk]])
                for vk in 1:length(valid_run_indices)]
            log_Z_mix_eff = _logsumexp(
                [log_alpha_eff[vk] + M * log1p(z0s[valid_run_indices[vk]])
                 for vk in 1:length(valid_run_indices)])
            log_pi_mix_eff = Vector{Vector{Float64}}(undef, length(valid_run_indices))
            for (vk, k) in enumerate(valid_run_indices)
                Ns_k = Ns_per_run[k]
                v = Vector{Float64}(undef, length(Ns_k))
                for j in 1:length(Ns_k)
                    v[j] = _log_mixture_prior_density(Ns_k[j],
                                                       log_normed_weights_eff,
                                                       log_z0s_valid)
                end
                log_pi_mix_eff[vk] = v
            end
        end

        # First pass over valid samples: find max log w̃ for stability.
        max_log = -Inf
        for (vk, k) in enumerate(valid_run_indices)
            log_α_k = log_alpha_eff[vk]
            Ns_k = Ns_per_run[k]
            Us_k = Us_per_run[k]
            log_ω_k = log_ω_per_run[k]
            log_pi_k = log_pi_mix_eff[vk]
            for j in 1:length(Ns_k)
                log_w_j = log_α_k + log_ω_k[j] + Ns_k[j] * (β * μ_target) -
                          β * Us_k[j] - log_pi_k[j]
                if log_w_j > max_log
                    max_log = log_w_j
                end
            end
        end

        if !isfinite(max_log)
            mean_E[i_μ, j_T]   = NaN
            mean_N[i_μ, j_T]   = NaN
            Cv_mu[i_μ, j_T]    = NaN
            ln_Z_rel[i_μ, j_T] = NaN
            ln_Z_abs[i_μ, j_T] = NaN
            N_eff_pooled[i_μ, j_T] = NaN
            continue
        end

        # Second pass: accumulate shifted moments.
        z = 0.0
        z2_shifted = 0.0
        e_sum = 0.0
        e2_sum = 0.0
        n_sum = 0.0
        en_sum = 0.0
        for (vk, k) in enumerate(valid_run_indices)
            log_α_k = log_alpha_eff[vk]
            Ns_k = Ns_per_run[k]
            Us_k = Us_per_run[k]
            log_ω_k = log_ω_per_run[k]
            log_pi_k = log_pi_mix_eff[vk]
            for j in 1:length(Ns_k)
                log_w_j = log_α_k + log_ω_k[j] + Ns_k[j] * (β * μ_target) -
                          β * Us_k[j] - log_pi_k[j]
                shifted = log_w_j - max_log
                w = exp(shifted)
                z += w
                z2_shifted += exp(2 * shifted)
                Uj = Us_k[j]
                Nj = Ns_k[j]
                e_sum += w * Uj
                e2_sum += w * Uj * Uj
                n_sum += w * Nj
                en_sum += w * Uj * Nj
            end
        end

        if z == 0.0 || !isfinite(z)
            mean_E[i_μ, j_T]   = NaN
            mean_N[i_μ, j_T]   = NaN
            Cv_mu[i_μ, j_T]    = NaN
            ln_Z_rel[i_μ, j_T] = NaN
            ln_Z_abs[i_μ, j_T] = NaN
            N_eff_pooled[i_μ, j_T] = NaN
            continue
        end

        e_avg = e_sum / z
        e2_avg = e2_sum / z
        n_avg = n_sum / z
        en_avg = en_sum / z
        var_e = e2_avg - e_avg^2
        cov_en = en_avg - e_avg * n_avg

        mean_E[i_μ, j_T]   = e_avg
        mean_N[i_μ, j_T]   = n_avg
        Cv_mu[i_μ, j_T]    = kb * β^2 * (var_e - μ_target * cov_en)
        # With the normalized-mixture denominator, LSE is ln Ξ_absolute directly.
        ln_Z_abs[i_μ, j_T] = max_log + log(z)
        # Bookkeeping: ln Z_relative reduces to single-run convention at K=1.
        ln_Z_rel[i_μ, j_T] = ln_Z_abs[i_μ, j_T] - log_Z_mix_eff
        N_eff_pooled[i_μ, j_T] = (z * z) / max(z2_shifted, eps(Float64))
    end

    base = (mean_E=mean_E,
            mean_N=mean_N,
            Cv_mu=Cv_mu,
            ln_Z_relative=ln_Z_rel,
            ln_Z_absolute=ln_Z_abs,
            N_eff_pooled=N_eff_pooled,
            N_eff_per_run=per_run_neff)
    return return_per_run ? merge(base, (per_run=per_run_stats,)) : base
end


end # module AnalysisTools