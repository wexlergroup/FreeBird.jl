"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV, Arrow
using Unitful

export read_output
export ωᵢ, partition_function, internal_energy, cv
export gc_thermodynamic_stats
export gc_thermodynamic_stats_fixed_N

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

After a finite number of NS iterations `n_iters` the recorded weights `ωᵢ` carry
only `1 − (K/(K+n_cull))^{n_iters}` of the prior volume (times `ω0`); the remainder
sits in the `K` surviving live walkers. Supplying `live_emax` (one vector of K live
walker energies per `N`) adds the live-set tail to each per-N evidence: each live
walker contributes weight `ω0 · (K/(K+n_cull))^{n_iters} / K` at its current energy.
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
- `V::typeof(1.0u"Å^3")`: accessible configurational volume (explicit; for
  surface systems this is *not* the simulation-box volume).
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
A `NamedTuple` `(Xi, mean_N, var_N, mean_U)`. Each field is a `Matrix{Float64}`
of size `(length(μ_grid), length(T_grid))` indexed `[i_μ, i_T]`. `Xi` is the
absolute grand partition function, `mean_N` is `⟨N⟩`, `var_N` is `⟨N²⟩ − ⟨N⟩²`,
and `mean_U` is `⟨E⟩` (grand-canonical, in eV).
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

    for (i, df) in enumerate(ns_outputs)
        # The N = 0 sector is the empty configuration: Z_NS^{(0)} = 1 by
        # definition. The corresponding DataFrame contents are ignored.
        if N_int[i] == 0
            log_Z_NS[i, :] .= 0.0
            mean_E_N[i, :] .= 0.0
            continue
        end

        ωi = ωᵢ(df.iter, n_walkers; n_cull=n_cull, ω0=ω0)
        Es = collect(Float64, df.emax)
        n_iters = length(df.iter)
        # Each live walker carries weight ω0 · (K/(K+n_cull))^n_iters / K,
        # accounting for the prior volume that finite termination leaves in
        # the live set.
        log_tail = log(ω0) + n_iters * log(n_walkers / (n_walkers + n_cull)) -
                   log(n_walkers)
        for (j, T) in enumerate(T_grid)
            β = 1.0 / (kb * ustrip(u"K", T))
            if live_emax === nothing
                log_terms = log.(ωi) .- β .* Es
                Es_all = Es
            else
                Es_live = collect(Float64, live_emax[i])
                log_terms = vcat(log.(ωi) .- β .* Es,
                                 log_tail .- β .* Es_live)
                Es_all = vcat(Es, Es_live)
            end
            max_log = maximum(log_terms)
            ws = exp.(log_terms .- max_log)
            sum_w = sum(ws)
            log_Z_NS[i, j] = max_log + log(sum_w)
            mean_E_N[i, j] = sum(ws .* Es_all) / sum_w
        end
    end

    log_fact = [_log_factorial(N) for N in N_int]

    Xi = Matrix{Float64}(undef, n_mu, n_T)
    mean_N = Matrix{Float64}(undef, n_mu, n_T)
    var_N = Matrix{Float64}(undef, n_mu, n_T)
    mean_U = Matrix{Float64}(undef, n_mu, n_T)

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
        end
    end

    return (Xi=Xi, mean_N=mean_N, var_N=var_N, mean_U=mean_U)
end


end # module AnalysisTools