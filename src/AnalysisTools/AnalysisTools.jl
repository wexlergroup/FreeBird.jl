"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV, Arrow
using Unitful

export read_output
export ωᵢ, log_ωᵢ, partition_function, internal_energy, cv
export gc_thermodynamic_stats
export gc_thermodynamic_stats_fixed_N
export microcanonical_entropy, caloric_derivatives, inflection_transitions
export transition_convergence
export gc_thermodynamic_stats_ideal_ref

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
\\omega_i = \\frac{C}{K+C} \\left(\\frac{K}{K+C}\\right)^{i-1}
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
    ωi = ω0 * (n_cull/(n_walkers+n_cull)) *
          (n_walkers/(n_walkers+n_cull)).^(iters .- 1)
    return ωi
end

"""
    log_ωᵢ(iters::Vector{Int}, n_walkers::Int; n_cull::Int=1, ω0::Float64=1.0)

Log of [`ωᵢ`](@ref), built directly in log space:

```math
\\log \\omega_i = \\log \\omega_0 + \\log\\frac{C}{K+C} + (i-1) \\log\\frac{K}{K+C}
```

Use this rather than `log.(ωᵢ(...))` anywhere the weights feed a log-sum-exp.
`ωᵢ` is a product of a factor slightly below one raised to the iteration
number, so it underflows to exactly `0.0` for `i ≳ 745·K/n_cull` — and
`log(0.0)` is `-Inf`, which silently drops those samples from the sum. They are
the *deepest*, lowest-energy samples, so what is lost is precisely the part that
dominates at low temperature. The scale is not exotic: at `K = 100` walkers with
`n_cull = 1` it starts at about 74,500 iterations.

This function is exact for any iteration count, and agrees with `log.(ωᵢ(...))`
everywhere `ωᵢ` has not underflowed.

# Arguments
- `iters::Vector{Int}`: The iteration numbers.
- `n_walkers::Int`: The number of walkers, \$K\$.
- `n_cull::Int`: The number of culled walkers, \$C\$. Default is 1.
- `ω0::Float64`: The initial \$\\omega\$ factor. Default is 1.0.

# Returns
- A vector of \$\\log \\omega\$ factors.
"""
function log_ωᵢ(iters::AbstractVector{Int}, n_walkers::Int; n_cull::Int=1, ω0::Float64=1.0)
    return (log(ω0) + log(n_cull / (n_walkers + n_cull))) .+
           (iters .- 1) .* log(n_walkers / (n_walkers + n_cull))
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

A `NamedTuple`. Its first three fields are `mean_E`, `cv`, `mean_N` in that
order, so `a, b, c = gc_thermodynamic_stats(...)` keeps working.

- `mean_E`: ⟨E⟩, the mean **bare** energy.
- `cv`: **`C_E = k_B β² [Var(E) − μ Cov(E,N)]`** — the thermodynamic heat
  capacity `(∂U/∂T)` at fixed μ and V. This is *the* heat capacity, and the
  default.
- `mean_N`: ⟨N⟩.
- `c_omega`: **`C_Ω = k_B β² Var(Ω)`**, `Ω = E − μN` — the fluctuation of the
  Hamiltonian actually sampled. Equal to `C_E` only at μ = 0; in general
  `C_Ω − C_E = −μ(∂⟨N⟩/∂T)_μ`, a difference comparable to the peak height
  itself near an order–disorder transition of a small adlayer. Reported rather
  than dropped because its peaks do locate transitions — but it is not
  `∂U/∂T`, and conflating the two is the reason both are named explicitly here.
  (For completeness: `C_Ω = T(∂S/∂T)_{μ,V}`, so it *is* a heat capacity of the
  open system in the `δQ = T dS` sense. Reporting `C_E` as "the" heat capacity
  is a choice of which response function to privilege, not a uniqueness claim.)
- `c_N`: **`C_N = k_B β² [Var(E) − Cov(E,N)²/Var(N)]`** — the part of the
  energy fluctuation uncorrelated with particle number. This is the standard
  grand-canonical-to-canonical ensemble conversion: applying
  `(∂U/∂T)_{N,V} = (∂U/∂T)_{μ,V} − (∂U/∂μ)_{T,V}(∂N/∂T)_{μ,V}/(∂N/∂μ)_{T,V}`
  to the four GC fluctuation identities cancels every μ-term and leaves exactly
  this expression, so `C_N` **is** the fixed-N heat capacity `C_{V,N}` measured
  in a variable-N ensemble (cf. Hill, *Statistical Mechanics*; Allen & Tildesley
  give the same formula for GCMC). It is the residual variance of E after linear
  regression on N, hence non-negative by Cauchy–Schwarz. Degenerates to
  `k_B β² Var(E)` when N does not fluctuate, where the projection is undefined
  rather than zero.
- `var_N`: Var(N).
"""
function gc_thermodynamic_stats(β::Float64,
                                 ωi::Vector{Float64},
                                 grand_energies::Vector{Float64},
                                 energies::Vector{Float64},
                                 numbers::Vector{Int},
                                 μ::Float64;
                                 kb::Float64=8.617333262e-5)
    # Public signature preserved. Callers holding linear weights get exactly the
    # behaviour they always did, underflow included — there is nothing to
    # recover once a weight has already reached 0.0. The DataFrame methods build
    # their weights with `log_ωᵢ` instead and call `_gc_stats_logw` directly.
    return _gc_stats_logw(β, log.(ωi), grand_energies, energies, numbers, μ; kb=kb)
end

"""
    _gc_stats_logw(β, log_ωi, grand_energies, energies, numbers, μ; kb)

Implementation of [`gc_thermodynamic_stats`](@ref) taking **log** weights, so the
weights can be constructed in log space by the caller and never round-trip
through a number that can underflow. Returns `(⟨E⟩, C_{V,μ}, ⟨N⟩)`.
"""
function _gc_stats_logw(β::Float64,
                        log_ωi::Vector{Float64},
                        grand_energies::Vector{Float64},
                        energies::Vector{Float64},
                        numbers::Vector{Int},
                        μ::Float64;
                        kb::Float64=8.617333262e-5)
    n = length(log_ωi)
    if n != length(grand_energies) || n != length(energies) || n != length(numbers)
        throw(DimensionMismatch("All input vectors must have the same length"))
    end
    if n == 0
        return _gc_stats_nan()
    end

    # Log-sum-exp for numerical stability
    log_terms = [log_ωi[i] - β * grand_energies[i] for i in 1:n]
    max_log = maximum(log_terms)

    z = 0.0
    u = 0.0   # ⟨E⟩
    u2 = 0.0  # ⟨E²⟩
    n_sum = 0.0  # ⟨N⟩
    n2_sum = 0.0 # ⟨N²⟩
    en_sum = 0.0 # ⟨EN⟩

    for i in 1:n
        w = exp(log_terms[i] - max_log)
        z += w
        u += w * energies[i]
        u2 += w * energies[i]^2
        n_sum += w * numbers[i]
        n2_sum += w * numbers[i]^2
        en_sum += w * energies[i] * numbers[i]
    end

    if z == 0.0
        return _gc_stats_nan()
    end

    u /= z
    u2 /= z
    n_avg = n_sum / z
    n2_avg = n2_sum / z
    en_avg = en_sum / z

    var_e = u2 - u^2
    var_n = n2_avg - n_avg^2
    cov_en = en_avg - u * n_avg

    cv, c_omega, c_N = _gc_heat_capacities(var_e, cov_en, var_n, n_avg, β, μ, kb)

    return (mean_E=u, cv=cv, mean_N=n_avg, c_omega=c_omega, c_N=c_N, var_N=var_n)
end

"""
    _gc_heat_capacities(var_e, cov_en, var_n, n_avg, β, μ, kb) -> (cv, c_omega, c_N)

The three grand-canonical heat-capacity definitions, in one place.

Each estimator in this module accumulates its weighted moments differently —
sequentially over samples, or vectorised over a (μ, T) grid — but they must
agree on what the *definitions* are. Keeping the formulas here is what stops
them drifting apart, which is precisely how the two heat capacities came to
disagree between the package and the prototype scripts in the first place.

- `cv` = `C_E = k_B β² [Var(E) − μ Cov(E,N)]`, the thermodynamic `(∂U/∂T)_{μ,V}`.
- `c_omega` = `C_Ω = k_B β² Var(Ω)` with `Ω = E − μN`, expanded as
  `Var(E) − 2μ Cov(E,N) + μ² Var(N)`. Equals `C_E` only at μ = 0.
- `c_N` = `k_B β² [Var(E) − Cov(E,N)²/Var(N)]`, the energy fluctuation with the
  part correlated with N projected out. When N does not fluctuate the
  projection is undefined rather than zero, so it degenerates to `k_B β² Var(E)`
  rather than dividing by ~0 and returning noise.
"""
function _gc_heat_capacities(var_e::Float64, cov_en::Float64, var_n::Float64,
                             n_avg::Float64, β::Float64, μ::Float64, kb::Float64)
    cv = kb * β^2 * (var_e - μ * cov_en)
    c_omega = kb * β^2 * (var_e - 2μ * cov_en + μ^2 * var_n)
    n_scale = max(1.0, abs(n_avg))
    c_N = var_n > 1e-12 * n_scale^2 ?
          kb * β^2 * (var_e - cov_en^2 / var_n) :
          kb * β^2 * var_e
    return cv, c_omega, c_N
end

"The all-NaN result, with the same fields as a successful one."
_gc_stats_nan() = (mean_E=NaN, cv=NaN, mean_N=NaN, c_omega=NaN, c_N=NaN, var_N=NaN)

"""
    _warn_if_reweighting(df, Es, Ns, μ)

Warn once if `μ` is not the chemical potential the run was sampled at.

Nothing here is an error: reweighting to a different μ is a supported use of a
GC-NS run, and is why `(Ω, E, N)` are all recorded rather than Ω alone. But the
estimator is only efficient near the μ the Ω-ladder was built at, and the
degradation is silent — the numbers come back looking fine — so it is worth
saying out loud.

The run's μ is recoverable from the data: `Ω = E − μN` exactly, for every
recorded row, so any row with `N ≠ 0` inverts to `μ_run = (E − Ω)/N`.
"""
function _warn_if_reweighting(df::DataFrame,
                              Es::Vector{Float64},
                              Ns::Vector{Int},
                              μ::Float64)
    (isempty(Es) || !hasproperty(df, :omega)) && return nothing
    ω_recorded = collect(Float64, df.omega)
    length(ω_recorded) == length(Es) || return nothing

    resid = maximum(abs.(ω_recorded .- (Es .- μ .* Ns)); init=0.0)
    tol = 1e-6 * max(1.0, maximum(abs, Es; init=1.0))
    resid <= tol && return nothing

    k = findfirst(!=(0), Ns)
    μ_run = k === nothing ? nothing : (Es[k] - ω_recorded[k]) / Ns[k]
    @warn "gc_thermodynamic_stats: μ is not the μ this run was sampled at. The result " *
          "is a post-hoc reweighting, which is supported but loses statistical " *
          "efficiency as |μ − μ_run| grows." requested_μ=μ μ_run=μ_run maxlog=1
    return nothing
end

"""
    gc_thermodynamic_stats(df::DataFrame, βs::Vector{Float64},
                           n_walkers::Int, μ::Float64;
                           n_cull::Int=1, ω0::Float64=1.0,
                           kb::Float64=8.617333262e-5)

Compute grand-canonical thermodynamic stats from a GC-NS output DataFrame.

The DataFrame must have columns `:iter`, `:omega`, `:energy`, `:num_particles`.
New sampler output also carries `energy_convention = "bare_E_v1"`; this makes
on-disk tables self-describing while preserving support for programmatically
constructed DataFrames.

Pass the surviving live set as `live_energies` / `live_numbers` to add the
live-walker contribution to the end of the recorded samples. Nested sampling
terminates after a finite number of iterations, and the prior volume left
unexplored at that point,

```math
X_f = \\left(\\frac{K}{K+C}\\right)^{i_f}
```

is carried entirely by the `K` surviving walkers. Omitting them truncates the
normalization at the deepest recorded sample, which biases ⟨E⟩, ⟨N⟩ and `Cv`
low-temperature-first — exactly where the live set dominates — and the shorter
the run, the worse it is. **Without these arguments the result is the truncated
estimate**, which is what this method returned unconditionally before, despite
its docstring.

# Arguments
- `df::DataFrame`: GC-NS output with columns `[:iter, :omega, :energy, :num_particles]`.
- `βs::Vector{Float64}`: Inverse temperatures at which to evaluate.
- `n_walkers::Int`: Number of walkers used in the NS run, \$K\$.
- `μ::Float64`: Chemical potential at which to evaluate. Ω is derived from
  `energy` and `num_particles` at this μ, so passing a μ other than the run's
  performs a post-hoc reweighting — supported, and warned about once, because
  its statistical efficiency degrades as `|μ − μ_run|` grows.
- `n_cull::Int=1`: Number of walkers culled per iteration, \$C\$.
- `ω0::Float64=1.0`: Initial phase-space volume.
- `live_energies=nothing`: Bare energies `E` of the surviving walkers, in the
  units of `df.energy`. Named `live_emax` in
  [`gc_thermodynamic_stats_ideal_ref`](@ref) and
  [`gc_thermodynamic_stats_fixed_N`](@ref), whose DataFrames call that column
  `emax`; here the energy column is `energy` and `omega` holds Ω.
- `live_numbers=nothing`: Particle numbers `N` of the surviving walkers. Their
  Ω is derived as `E - μN`, so the live set does not need to have been recorded
  at this μ.
- `kb::Float64`: Boltzmann constant (default: eV/K).

# Returns

A `NamedTuple` of vectors, one entry per β. Its first three fields are `mean_E`, `cv`, `mean_N` in that
order, so `a, b, c = gc_thermodynamic_stats(...)` keeps working.

- `mean_E`: ⟨E⟩ at each β, the mean **bare** energy.
- `cv`: **`C_E = k_B β² [Var(E) − μ Cov(E,N)]`** — the thermodynamic heat
  capacity `(∂U/∂T)` at fixed μ and V. This is *the* heat capacity, and the
  default.
- `mean_N`: ⟨N⟩.
- `c_omega`: **`C_Ω = k_B β² Var(Ω)`**, `Ω = E − μN` — the fluctuation of the
  Hamiltonian actually sampled. Equal to `C_E` only at μ = 0; in general
  `C_Ω − C_E = −μ(∂⟨N⟩/∂T)_μ`, a difference comparable to the peak height
  itself near an order–disorder transition of a small adlayer. Reported rather
  than dropped because its peaks do locate transitions — but it is not
  `∂U/∂T`, and conflating the two is the reason both are named explicitly here.
  (For completeness: `C_Ω = T(∂S/∂T)_{μ,V}`, so it *is* a heat capacity of the
  open system in the `δQ = T dS` sense. Reporting `C_E` as "the" heat capacity
  is a choice of which response function to privilege, not a uniqueness claim.)
- `c_N`: **`C_N = k_B β² [Var(E) − Cov(E,N)²/Var(N)]`** — the part of the
  energy fluctuation uncorrelated with particle number. This is the standard
  grand-canonical-to-canonical ensemble conversion: applying
  `(∂U/∂T)_{N,V} = (∂U/∂T)_{μ,V} − (∂U/∂μ)_{T,V}(∂N/∂T)_{μ,V}/(∂N/∂μ)_{T,V}`
  to the four GC fluctuation identities cancels every μ-term and leaves exactly
  this expression, so `C_N` **is** the fixed-N heat capacity `C_{V,N}` measured
  in a variable-N ensemble (cf. Hill, *Statistical Mechanics*; Allen & Tildesley
  give the same formula for GCMC). It is the residual variance of E after linear
  regression on N, hence non-negative by Cauchy–Schwarz. Degenerates to
  `k_B β² Var(E)` when N does not fluctuate, where the projection is undefined
  rather than zero.
- `var_N`: Var(N).
"""
function gc_thermodynamic_stats(df::DataFrame,
                                 βs::Vector{Float64},
                                 n_walkers::Int,
                                 μ::Float64;
                                 n_cull::Int=1,
                                 ω0::Float64=1.0,
                                 live_energies=nothing,
                                 live_numbers=nothing,
                                 kb::Float64=8.617333262e-5)
    log_ωi = log_ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
    Es = collect(Float64, df.energy)
    Ns = collect(Int, df.num_particles)

    # Ω is *derived* at the requested μ rather than read from `df.omega`. At the
    # run's own μ these are the same numbers — the sampler records
    # `omega = energy - mu * num_particles` from exactly these columns — so the
    # ordinary call is unchanged. It matters when μ differs from the run's:
    # that is the post-hoc reweighting the recorded (Ω, E, N) triple exists to
    # make possible, and reading `df.omega` there would weight the samples on
    # the run's Ω-ladder while correcting `Cv` at the requested μ, mixing two
    # ensembles in one estimate.
    grand_es = Es .- μ .* Ns
    _warn_if_reweighting(df, Es, Ns, μ)

    if live_energies !== nothing && !isempty(live_energies)
        if live_numbers === nothing
            throw(ArgumentError("live_energies given without live_numbers"))
        end
        if length(live_energies) != length(live_numbers)
            throw(DimensionMismatch("live_energies and live_numbers must have the same length"))
        end
        # Residual prior volume after the last recorded iteration, split
        # uniformly over the K surviving walkers. No ω0 factor here: ω0 only
        # rescales the dead-sample shell weights. With the normalized default
        # ω0 = 1, the dead weights sum to 1 − X_f and the tail closes Σw = 1.
        # Same construction as gc_thermodynamic_stats_ideal_ref.
        n_iters = isempty(df.iter) ? 0 : maximum(df.iter)
        log_tail = n_iters * log(n_walkers / (n_walkers + n_cull)) - log(n_walkers)

        Es_live = collect(Float64, live_energies)
        Ns_live = collect(Int, live_numbers)

        log_ωi = vcat(log_ωi, fill(log_tail, length(Es_live)))
        grand_es = vcat(grand_es, Es_live .- μ .* Ns_live)
        Es = vcat(Es, Es_live)
        Ns = vcat(Ns, Ns_live)
    end

    nβ = length(βs)
    mean_Es = Vector{Float64}(undef, nβ)
    Cvs = Vector{Float64}(undef, nβ)
    mean_Ns = Vector{Float64}(undef, nβ)
    c_omegas = Vector{Float64}(undef, nβ)
    c_Ns = Vector{Float64}(undef, nβ)
    var_Ns = Vector{Float64}(undef, nβ)

    Threads.@threads for (i, b) in collect(enumerate(βs))
        r = _gc_stats_logw(b, log_ωi, grand_es, Es, Ns, μ; kb=kb)
        mean_Es[i] = r.mean_E
        Cvs[i] = r.cv
        mean_Ns[i] = r.mean_N
        c_omegas[i] = r.c_omega
        c_Ns[i] = r.c_N
        var_Ns[i] = r.var_N
    end

    return (mean_E=mean_Es, cv=Cvs, mean_N=mean_Ns,
            c_omega=c_omegas, c_N=c_Ns, var_N=var_Ns)
end

"""
    gc_thermodynamic_stats_ideal_ref(df::DataFrame, n_sites::Int, z0::Float64,
                                     μs::Vector{Float64}, Ts::Vector{Float64},
                                     n_walkers::Int;
                                     n_cull::Int=1, ω0::Float64=1.0,
                                     live_emax::Union{Nothing,Vector{Float64}}=nothing,
                                     live_numbers::Union{Nothing,Vector{Int}}=nothing,
                                     kb::Float64=8.617333262e-5)

Assemble grand-canonical thermodynamics on a (μ, T) grid from a single
ideal-gas-referenced nested sampling run (`ideal_gas_referenced_nested_sampling`
in `SamplingSchemes`).

The run samples the ideal-lattice-gas prior at reference fugacity `z0`
(configuration weight `z0^N`, total prior mass `(1 + z0)^M` on `M = n_sites`
sites) and culls by energy alone. The absolute grand partition function at any
target chemical potential μ and temperature T is then

```math
\\Xi(\\mu, T) = (1 + z_0)^M \\sum_j \\omega_j \\left(\\frac{z}{z_0}\\right)^{N_j}
               e^{-\\beta E_j}, \\qquad z = e^{\\beta\\mu},
```

where the sum runs over culled walkers (plus, when `live_emax`/`live_numbers`
are given, the surviving live walkers, each with residual weight
`(K/(K+n_cull))^{n_iters} / K` — no `ω0` factor, since `ω0` rescales the
dead-sample shell weights only). There is no thermal-wavelength factor: a
lattice gas has no momentum degrees of freedom, so `z = exp(βμ)` directly.
All sums are evaluated with the log-sum-exp trick, and Ξ is returned as
`log Ξ` to avoid overflow at low temperature.

For a fully normalized absolute Ξ, use the default `ω0 = 1.0` together with the
live-walker tail. The one-based Skilling shell weights then sum to
`1 − (K/(K+n_cull))^{n_iters}` and the tail supplies the remainder, so `Σω = 1`
exactly. Omitting the live set neglects the residual prior volume. Ratio
observables (`mean_N`, `var_N`, `mean_U`) are insensitive to a common `ω0`
when no live tail is supplied.

The reweighting factor `(z/z0)^{N_j}` is pure importance sampling in μ: its
reliability at each grid point is reported by the Kish effective sample size
`N_eff = (Σ w)² / Σ w²`, which collapses as `|βμ − ln z0|` grows beyond
roughly `1/√Var(N)`. Treat grid points with small `N_eff` (≲ 100) as
unreliable and re-run with z0 closer to the target fugacity.

# Arguments
- `df::DataFrame`: NS output with columns `[:iter, :emax, :num_particles]`.
- `n_sites::Int`: Number of lattice sites M.
- `z0::Float64`: Reference fugacity of the prior used in the run (must match!).
- `μs::Vector{Float64}`: Chemical potential grid (same energy units as `df.emax`, e.g. eV).
- `Ts::Vector{Float64}`: Temperature grid in K.
- `n_walkers::Int`: Number of walkers K used in the NS run.
- `n_cull::Int=1`: Number of walkers culled per iteration.
- `ω0::Float64=1.0`: Initial phase-space volume factor.
- `live_emax::Union{Nothing,Vector{Float64}}=nothing`: Energies of the surviving live walkers.
- `live_numbers::Union{Nothing,Vector{Int}}=nothing`: Particle counts of the surviving live walkers.
- `kb::Float64`: Boltzmann constant (default: eV/K).

# Returns
A `NamedTuple` of `Matrix{Float64}` of size `(length(μs), length(Ts))`,
indexed `[i_μ, i_T]`:
- `logXi`: Natural log of the absolute grand partition function.
- `mean_N`: Mean particle number ⟨N⟩.
- `var_N`: Particle-number variance ⟨N²⟩ − ⟨N⟩².
- `mean_U`: Mean configurational energy ⟨E⟩.
- `cv`: `C_E = k_B β² [Var(E) − μ Cov(E,N)]`, the thermodynamic heat capacity
  `(∂U/∂T)_{μ,V}` — the default heat capacity, and the one to quote.
- `c_omega`: `C_Ω = k_B β² Var(E − μN)`, the fluctuation of the Hamiltonian a
  fixed-μ run would have sampled. Equals `C_E` only at μ = 0.
- `c_N`: `k_B β² [Var(E) − Cov(E,N)²/Var(N)]`, the energy fluctuation with the
  part correlated with particle number projected out.

All three come from `_gc_heat_capacities`, shared with `gc_thermodynamic_stats`
so the definitions cannot drift between the two estimators. They inherit the
reweighting caveat above: at grid points where `N_eff` has collapsed they are
as unreliable as everything else there, and more so, being second moments.
- `N_eff`: Kish effective sample size of the reweighted estimate.
"""
function gc_thermodynamic_stats_ideal_ref(df::DataFrame,
                                          n_sites::Int,
                                          z0::Float64,
                                          μs::Vector{Float64},
                                          Ts::Vector{Float64},
                                          n_walkers::Int;
                                          n_cull::Int=1,
                                          ω0::Float64=1.0,
                                          live_emax::Union{Nothing,Vector{Float64}}=nothing,
                                          live_numbers::Union{Nothing,Vector{Int}}=nothing,
                                          kb::Float64=8.617333262e-5)
    if z0 <= 0.0
        throw(ArgumentError("z0 must be positive"))
    end
    if n_sites <= 0
        throw(ArgumentError("n_sites must be positive"))
    end
    if (live_emax === nothing) != (live_numbers === nothing)
        throw(ArgumentError("live_emax and live_numbers must be provided together"))
    end
    if live_emax !== nothing && length(live_emax) != length(live_numbers)
        throw(DimensionMismatch("live_emax and live_numbers must have the same length"))
    end
    n_dead = nrow(df)
    if n_dead == 0 && (live_emax === nothing || isempty(live_emax))
        throw(ArgumentError("df is empty and no live walkers were provided"))
    end

    # β- and μ-independent per-sample log prior-volume weights, built directly
    # in log space: the linear ωᵢ underflows to 0.0 once iter ≳ 745·K/n_cull,
    # which silently zeroes the deepest (lowest-energy) samples on large lattices
    log_w0 = n_dead > 0 ?
        (log(ω0) + log(n_cull / (n_walkers + n_cull))) .+
        (Vector{Float64}(df.iter) .- 1) .* log(n_walkers / (n_walkers + n_cull)) : Float64[]
    Es = n_dead > 0 ? Vector{Float64}(df.emax) : Float64[]
    Ns = n_dead > 0 ? Vector{Float64}(df.num_particles) : Float64[]

    if live_emax !== nothing && !isempty(live_emax)
        # Residual prior volume after the last recorded iteration, split
        # uniformly over the K surviving walkers. No ω0 factor here: ω0
        # rescales the dead-sample shell weights, not the residual volume
        # X_n = (K/(K+n_cull))^n. With ω0 = 1 the dead weights sum to
        # 1 − X_n and the tail closes the identity Σw = 1 exactly.
        n_iters = n_dead > 0 ? maximum(df.iter) : 0
        log_tail = n_iters * log(n_walkers / (n_walkers + n_cull)) - log(n_walkers)
        log_w0 = vcat(log_w0, fill(log_tail, length(live_emax)))
        Es = vcat(Es, live_emax)
        Ns = vcat(Ns, Float64.(live_numbers))
    end

    log_prior_mass = n_sites * log1p(z0)
    log_z0 = log(z0)

    n_mu = length(μs)
    n_T = length(Ts)
    logXi = Matrix{Float64}(undef, n_mu, n_T)
    mean_N = Matrix{Float64}(undef, n_mu, n_T)
    var_N = Matrix{Float64}(undef, n_mu, n_T)
    mean_U = Matrix{Float64}(undef, n_mu, n_T)
    N_eff = Matrix{Float64}(undef, n_mu, n_T)
    cv = Matrix{Float64}(undef, n_mu, n_T)
    c_omega = Matrix{Float64}(undef, n_mu, n_T)
    c_N = Matrix{Float64}(undef, n_mu, n_T)

    Threads.@threads for j in 1:n_T
        β = 1.0 / (kb * Ts[j])
        for i in 1:n_mu
            s = β * μs[i] - log_z0
            log_w = log_w0 .+ s .* Ns .- β .* Es
            max_log = maximum(log_w)
            ws = exp.(log_w .- max_log)
            sum_w = sum(ws)
            logXi[i, j] = log_prior_mass + max_log + log(sum_w)
            n_avg = sum(ws .* Ns) / sum_w
            u = sum(ws .* Es) / sum_w
            vN = sum(ws .* Ns .^ 2) / sum_w - n_avg^2
            vE = sum(ws .* Es .^ 2) / sum_w - u^2
            cEN = sum(ws .* Es .* Ns) / sum_w - u * n_avg
            mean_N[i, j] = n_avg
            var_N[i, j] = vN
            mean_U[i, j] = u
            N_eff[i, j] = sum_w^2 / sum(abs2, ws)
            cv[i, j], c_omega[i, j], c_N[i, j] =
                _gc_heat_capacities(vE, cEN, vN, n_avg, β, μs[i], kb)
        end
    end

    return (logXi=logXi, mean_N=mean_N, var_N=var_N, mean_U=mean_U, N_eff=N_eff,
            cv=cv, c_omega=c_omega, c_N=c_N)
end


"""
    _thermal_wavelength(atomic_mass::typeof(1.0u"u"), T::Unitful.Temperature) -> typeof(1.0u"Å")

Compute the thermal de Broglie wavelength `Λ = h / sqrt(2π m k_B T)` in Å.
Internal helper for `gc_thermodynamic_stats_fixed_N`.
"""
function _thermal_wavelength(atomic_mass::typeof(1.0u"u"), T::Unitful.Temperature)
    return uconvert(u"Å", Unitful.h / sqrt(2π * atomic_mass * Unitful.k * T))
end

"""
    _log_factorial(n::Integer) -> Float64

Return `log(n!)` for small non-negative `n`, computed by summing `log(k)` to
keep the result exact in floating point for `n` up to several dozen.
"""
function _log_factorial(n::Integer)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    n == 0 && return 0.0
    return sum(log, 1:n)
end

"""
    gc_thermodynamic_stats_fixed_N(ns_outputs, N_values, V, atomic_mass, μ_grid, T_grid;
                                   n_walkers=120, n_cull=1, ω0=1.0,
                                   live_emax=nothing, kb=8.617333262e-5)

Compute grand-canonical thermodynamic averages from a stack of canonical
nested-sampling outputs, one per fixed particle number `N`.

For each `N`, the canonical NS evidence at inverse temperature `β` is

```math
Z_{\\mathrm{NS}}^{(N)}(\\beta) = \\sum_i \\omega_i \\exp(-\\beta E_i)
```

— the prior-volume-normalized configurational integral. The absolute
configurational partition function is `Z_N^{config}(\\beta) = V^N \\cdot Z_{NS}^{(N)}(\\beta)`.
The grand partition function is assembled as

```math
\\Xi(\\mu, T) = \\sum_N \\frac{(zV)^N}{N!}\\, Z_{\\mathrm{NS}}^{(N)}(\\beta),
\\qquad z = \\frac{\\exp(\\beta\\mu)}{\\Lambda(T)^3}
```

with the thermal wavelength `Λ(T) = h / sqrt(2π m k_B T)` computed from `atomic_mass`.
The sum runs over the supplied `N_values`, which must include `0`. The `N=0` sector
is treated specially: `Z_{NS}^{(0)} = 1` by definition (the empty configuration has no
spatial integral), so the corresponding DataFrame contents are ignored. Truncation
error at the upper end of `N_values` is bounded by the tail of `(zV)^N / N!` for the
largest `⟨N⟩` requested.

A log-sum-exp pass is used both inside each per-N evidence and across the
grand sum for numerical stability.

## Live-set tail correction

After a finite number of NS iterations `n_iters` the recorded weights `ωᵢ` sum to
`ω0 · (1 − r^{n_iters})` with `r = K/(K+n_cull)`; the remaining prior volume,
`X_f = r^{n_iters}`, sits in the `K` surviving live walkers. Supplying `live_emax`
(one vector of K live walker energies per `N`) adds that tail to each per-N
evidence: each live walker contributes weight `X_f / K` at its current energy —
with **no** `ω0` factor, since `ω0` rescales the dead-sample shell weights and not
the residual volume. At the normalized default `ω0 = 1` the dead weights sum to
exactly `1 − X_f` and the tail closes `Σw = 1`.

When omitted, the live-set tail is neglected — for ratio observables (`⟨N⟩`, `⟨U⟩`)
the resulting bias is small but visible at low T or shallow NS; for the absolute
`Ξ` it appears as a uniform-in-N prefactor that does not cancel.

# Arguments
- `ns_outputs::AbstractVector{<:DataFrame}`: one canonical-NS output per `N`,
  each with columns `[:iter, :emax]` (matching the schema produced by
  `nested_sampling`). The entry corresponding to `N=0` is ignored and may be
  any DataFrame (e.g., `DataFrame(iter=Int[], emax=Float64[])`).
- `N_values::AbstractVector{<:Integer}`: particle counts corresponding to each
  DataFrame. Must include `0`; `length(N_values) == length(ns_outputs)`.
- `V::typeof(1.0u"Å^3")`: the simulation-box volume, i.e. the NS prior volume
  per particle (NS samples positions uniformly over the box). This is what
  closes `Z_N^{config} = V^N · Z_{NS}^{(N)}`, so it must *not* be reduced to an
  accessible or adsorption-region sub-volume. For surface systems the
  substrate's volume exclusion is captured by the Boltzmann factor inside
  `Z_{NS}^{(N)}`, not by shrinking `V`.
- `atomic_mass::typeof(1.0u"u")`: per-atom mass for `Λ(T)`.
- `μ_grid::AbstractVector{<:typeof(1.0u"eV")}`: chemical potentials.
- `T_grid::AbstractVector{<:Unitful.Temperature}`: temperatures.

# Keyword arguments
- `n_walkers::Int=120`: number of NS walkers used to produce each DataFrame.
  Must be uniform across `ns_outputs`.
- `n_cull::Int=1`: NS culls per iteration.
- `ω0::Float64=1.0`: initial prior weight, passed to `ωᵢ`.
- `live_emax::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing`:
  when supplied, one vector of `K = n_walkers` live walker energies (in eV) per
  `N`. The entry for `N=0` is ignored. See "Live-set tail correction" above.
- `kb::Float64`: Boltzmann constant in eV/K.

# Returns
A `NamedTuple` `(Xi, mean_N, var_N, mean_U, cv, c_omega, c_N)`. Each field is a
`Matrix{Float64}` of size `(length(μ_grid), length(T_grid))` indexed
`[i_μ, i_T]`. `Xi` is the absolute grand partition function, `mean_N` is `⟨N⟩`,
`var_N` is `⟨N²⟩ − ⟨N⟩²`, and `mean_U` is `⟨E⟩` (grand-canonical, in eV).

`cv` (`C_E`, the thermodynamic heat capacity and the one to quote), `c_omega`
(`C_Ω`) and `c_N` come from `_gc_heat_capacities`, shared with the other two
estimators so the definitions cannot drift. Their `Var(E)` is assembled across
the `N` sectors by the law of total variance — the per-sector `⟨E²⟩` is carried
forward alongside `⟨E⟩`, so the fluctuation *within* each fixed-`N` run is
included rather than only the spread of the sector means.
"""
function gc_thermodynamic_stats_fixed_N(
    ns_outputs::AbstractVector{<:DataFrame},
    N_values::AbstractVector{<:Integer},
    V::typeof(1.0u"Å^3"),
    atomic_mass::typeof(1.0u"u"),
    μ_grid::AbstractVector{<:typeof(1.0u"eV")},
    T_grid::AbstractVector{<:Unitful.Temperature};
    n_walkers::Int=120,
    n_cull::Int=1,
    ω0::Float64=1.0,
    live_emax::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing,
    kb::Float64=8.617333262e-5,
)
    if length(ns_outputs) != length(N_values)
        throw(DimensionMismatch("ns_outputs and N_values must have the same length"))
    end
    if !(0 in N_values)
        throw(ArgumentError("N_values must include 0 (the empty configuration)"))
    end
    if live_emax !== nothing && length(live_emax) != length(N_values)
        throw(DimensionMismatch("live_emax and N_values must have the same length"))
    end

    N_int = collect(Int, N_values)
    n_N = length(N_int)
    n_mu = length(μ_grid)
    n_T = length(T_grid)
    V_val = ustrip(u"Å^3", V)

    log_Z_NS = Matrix{Float64}(undef, n_N, n_T)
    mean_E_N = Matrix{Float64}(undef, n_N, n_T)
    # ⟨E²⟩ within each N sector. Needed for the grand-canonical Var(E), which by
    # the law of total variance is E_N[Var(E|N)] + Var_N[⟨E|N⟩] — the first term
    # is invisible if only the per-N means are carried forward.
    mean_E2_N = Matrix{Float64}(undef, n_N, n_T)

    for (i, df) in enumerate(ns_outputs)
        # The N = 0 sector is the empty configuration: Z_NS^{(0)} = 1 by
        # definition. The corresponding DataFrame contents are ignored.
        if N_int[i] == 0
            log_Z_NS[i, :] .= 0.0
            mean_E_N[i, :] .= 0.0
            mean_E2_N[i, :] .= 0.0
            continue
        end

        log_ωi = log_ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
        Es = collect(Float64, df.emax)
        # Residual prior volume after the last recorded iteration, split
        # uniformly over the K surviving walkers. No ω0 factor here: ω0
        # rescales the dead-sample shell weights, not the residual volume X_f.
        # With ω0 = 1 the dead weights sum to 1 − X_f and the tail closes
        # Σw = 1 exactly. Identical construction to gc_thermodynamic_stats and
        # gc_thermodynamic_stats_ideal_ref, which is what the docstrings above
        # claim all three share; this one used to carry a stray log(ω0), which
        # overweighted the tail by (K+C)/K, and indexed n_iters by `length`
        # rather than `maximum` of df.iter.
        n_iters = isempty(df.iter) ? 0 : maximum(df.iter)
        log_tail = n_iters * log(n_walkers / (n_walkers + n_cull)) - log(n_walkers)
        for (j, T) in enumerate(T_grid)
            β = 1.0 / (kb * ustrip(u"K", T))
            if live_emax === nothing
                log_terms = log_ωi .- β .* Es
                Es_all = Es
            else
                Es_live = collect(Float64, live_emax[i])
                log_terms = vcat(log_ωi .- β .* Es,
                                 log_tail .- β .* Es_live)
                Es_all = vcat(Es, Es_live)
            end
            max_log = maximum(log_terms)
            ws = exp.(log_terms .- max_log)
            sum_w = sum(ws)
            log_Z_NS[i, j] = max_log + log(sum_w)
            mean_E_N[i, j] = sum(ws .* Es_all) / sum_w
            mean_E2_N[i, j] = sum(ws .* Es_all .^ 2) / sum_w
        end
    end

    log_fact = [_log_factorial(N) for N in N_int]

    Xi = Matrix{Float64}(undef, n_mu, n_T)
    mean_N = Matrix{Float64}(undef, n_mu, n_T)
    var_N = Matrix{Float64}(undef, n_mu, n_T)
    mean_U = Matrix{Float64}(undef, n_mu, n_T)
    cv = Matrix{Float64}(undef, n_mu, n_T)
    c_omega = Matrix{Float64}(undef, n_mu, n_T)
    c_N = Matrix{Float64}(undef, n_mu, n_T)

    for (j, T) in enumerate(T_grid)
        β = 1.0 / (kb * ustrip(u"K", T))
        Λ_val = ustrip(u"Å", _thermal_wavelength(atomic_mass, T))
        log_V_over_Λ3 = log(V_val) - 3 * log(Λ_val)
        for (k, μ) in enumerate(μ_grid)
            μ_val = ustrip(u"eV", μ)
            log_zV = β * μ_val + log_V_over_Λ3

            log_w = log_Z_NS[:, j] .+ N_int .* log_zV .- log_fact
            max_log = maximum(log_w)
            ws = exp.(log_w .- max_log)
            sum_w = sum(ws)

            Xi[k, j] = exp(max_log) * sum_w
            mean_N[k, j] = sum(ws .* N_int) / sum_w
            mean_N2 = sum(ws .* (N_int .^ 2)) / sum_w
            var_N[k, j] = mean_N2 - mean_N[k, j]^2
            mean_U[k, j] = sum(ws .* view(mean_E_N, :, j)) / sum_w

            # Grand-canonical second moments assembled across the N sectors.
            # ⟨E²⟩ uses the per-sector ⟨E²⟩ rather than ⟨E⟩², which is what
            # keeps the within-sector fluctuation in the total.
            mean_E2 = sum(ws .* view(mean_E2_N, :, j)) / sum_w
            mean_EN = sum(ws .* N_int .* view(mean_E_N, :, j)) / sum_w
            var_e = mean_E2 - mean_U[k, j]^2
            cov_en = mean_EN - mean_U[k, j] * mean_N[k, j]
            cv[k, j], c_omega[k, j], c_N[k, j] =
                _gc_heat_capacities(var_e, cov_en, var_N[k, j], mean_N[k, j],
                                    β, μ_val, kb)
        end
    end

    return (Xi=Xi, mean_N=mean_N, var_N=var_N, mean_U=mean_U,
            cv=cv, c_omega=c_omega, c_N=c_N)
end


include("microcanonical_inflection.jl")

end # module AnalysisTools
