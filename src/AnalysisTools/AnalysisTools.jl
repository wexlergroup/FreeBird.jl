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
`(K/(K+n_cull))^{n_iters} / K` — no `ω0` factor, since `ω0` corrects the
dead-sample shell weights only). There is no thermal-wavelength factor: a
lattice gas has no momentum degrees of freedom, so `z = exp(βμ)` directly.
All sums are evaluated with the log-sum-exp trick, and Ξ is returned as
`log Ξ` to avoid overflow at low temperature.

For a fully normalized absolute Ξ, pass `ω0 = (n_walkers + n_cull)/n_walkers`
(Skilling weights) together with the live-walker tail — the dead weights then
sum to `1 − (K/(K+n_cull))^{n_iters}` and the tail supplies the remainder, so
`Σω = 1` exactly. The default `ω0 = 1.0` underestimates Ξ by a factor
`K/(K+n_cull)` and neglects the tail (see the `ωᵢ` conventions). Ratio
observables (`mean_N`, `var_N`, `mean_U`) are insensitive to `ω0`.

The reweighting factor `(z/z0)^{N_j}` is pure importance sampling in μ: its
reliability at each grid point is reported by the Kish effective sample size
`N_eff = (Σ w)² / Σ w²`, which collapses as `|βμ − ln z0|` grows beyond
roughly `1/√Var(N)`. Treat grid points with small `N_eff` (≲ 100) as
unreliable and re-run with z0 closer to the target fugacity.

Athermal (hard-core) models sampled with the finite-J recipe (see the
`GenericLatticeHamiltonian` docstring) are evaluated here at a temperature low
enough that `β·J` exceeds the ladder depth in nats by ~40 while `β·δ ≪ 1`
(δ = the run's `energy_perturbation`): the violating shells' `e^{-βE_j}` are
then exact zeros, the retained shells' deviate from 1 only by `O(β·δ)`, and
every observable is a function of `z = exp(βμ)` alone (`mean_U` reports the
residual tie-breaking noise, ~δ).

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
- `observable_cols::AbstractVector{Symbol}=Symbol[]`: Names of extra `df`
  columns (recorded per dead point, e.g. via the `observables` keyword of
  `ideal_gas_referenced_nested_sampling`) whose reweighted averages
  ⟨A⟩(μ, T) are returned in `observables`. Requesting `:num_particles` or
  `:emax` legitimately reproduces `mean_N`/`mean_U` (a self-consistency
  check).
- `live_observables=nothing`: Required whenever `observable_cols` is
  non-empty and a live-set tail is supplied — a
  `Dict{Symbol,<:AbstractVector{<:Real}}` with exactly the
  `observable_cols` keys, holding each observable's values on the surviving
  live walkers (one entry per `live_emax` entry, same order). Omitting the
  tail values would silently deflate every ⟨A⟩, so this is an error rather
  than a default.
- `kb::Float64`: Boltzmann constant (default: eV/K).

# Returns
A `NamedTuple`; every field except `observables` is a `Matrix{Float64}` of
size `(length(μs), length(Ts))`, indexed `[i_μ, i_T]`:
- `logXi`: Natural log of the absolute grand partition function.
- `mean_N`: Mean particle number ⟨N⟩.
- `var_N`: Particle-number variance ⟨N²⟩ − ⟨N⟩².
- `mean_U`: Mean configurational energy ⟨E⟩.
- `N_eff`: Kish effective sample size of the reweighted estimate.
- `var_U`: Energy variance ⟨E²⟩ − ⟨E⟩². The grand-canonical heat capacity
  follows as `C = k_B β² (var_U − μ · cov_UN)`, matching the Ω-sorted
  `gc_thermodynamic_stats` formula.
- `cov_UN`: Energy–particle-number covariance ⟨EN⟩ − ⟨E⟩⟨N⟩.
- `observables`: `Dict{Symbol,Matrix{Float64}}` mapping each requested
  column to its reweighted average ⟨A⟩(μ, T); empty when no columns are
  requested.
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
                                          observable_cols::AbstractVector{Symbol}=Symbol[],
                                          live_observables::Union{Nothing,AbstractDict{Symbol,<:AbstractVector{<:Real}}}=nothing,
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
    allunique(observable_cols) || throw(ArgumentError(
        "observable_cols must not contain duplicate entries"))
    for col in observable_cols
        hasproperty(df, col) || throw(ArgumentError(
            "observable_cols: df has no column :$col"))
        any(ismissing, df[!, col]) && throw(ArgumentError(
            "observable_cols: df column :$col contains missing values"))
    end
    if live_observables !== nothing && live_emax === nothing
        throw(ArgumentError(
            "live_observables requires live_emax/live_numbers (the live-set tail)"))
    end
    if live_observables !== nothing && isempty(observable_cols)
        throw(ArgumentError(
            "live_observables was given but observable_cols is empty"))
    end
    if !isempty(observable_cols) && live_emax !== nothing
        if live_observables === nothing
            throw(ArgumentError(
                "observable_cols with a live-set tail requires live_observables " *
                "(the per-column values of the surviving live walkers); omitting " *
                "them would silently deflate the tail's contribution to every " *
                "observable average"))
        end
        if Set(keys(live_observables)) != Set(observable_cols)
            throw(ArgumentError(
                "live_observables keys must match observable_cols exactly"))
        end
        for col in observable_cols
            if length(live_observables[col]) != length(live_emax)
                throw(DimensionMismatch(
                    "live_observables[:$col] must have one entry per live " *
                    "walker (length(live_emax) = $(length(live_emax)))"))
            end
        end
    end

    # β- and μ-independent per-sample log prior-volume weights, built directly
    # in log space: the linear ωᵢ underflows to 0.0 once iter ≳ 745·K/n_cull,
    # which silently zeroes the deepest (lowest-energy) samples on large lattices
    log_w0 = n_dead > 0 ?
        (log(ω0) + log(n_cull / (n_walkers + n_cull))) .+
        Vector{Float64}(df.iter) .* log(n_walkers / (n_walkers + n_cull)) : Float64[]
    Es = n_dead > 0 ? Vector{Float64}(df.emax) : Float64[]
    Ns = n_dead > 0 ? Vector{Float64}(df.num_particles) : Float64[]

    if live_emax !== nothing && !isempty(live_emax)
        # Residual prior volume after the last recorded iteration, split
        # uniformly over the K surviving walkers. No ω0 factor here: ω0
        # corrects the dead-sample shell weights, not the residual volume
        # X_n = (K/(K+n_cull))^n — with ω0 = (K+n_cull)/K the dead weights sum
        # to 1 − X_n and the tail closes the identity Σw = 1 exactly
        n_iters = n_dead > 0 ? maximum(df.iter) : 0
        log_tail = n_iters * log(n_walkers / (n_walkers + n_cull)) - log(n_walkers)
        log_w0 = vcat(log_w0, fill(log_tail, length(live_emax)))
        Es = vcat(Es, live_emax)
        Ns = vcat(Ns, Float64.(live_numbers))
    end

    # Per-column observable values aligned with (Es, Ns): dead points first,
    # then the live-set tail in the same order as live_emax
    As = Dict{Symbol,Vector{Float64}}()
    for col in observable_cols
        A = n_dead > 0 ? Vector{Float64}(df[!, col]) : Float64[]
        if live_emax !== nothing && !isempty(live_emax)
            A = vcat(A, Vector{Float64}(live_observables[col]))
        end
        As[col] = A
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
    var_U = Matrix{Float64}(undef, n_mu, n_T)
    cov_UN = Matrix{Float64}(undef, n_mu, n_T)
    obs_out = Dict{Symbol,Matrix{Float64}}(
        col => Matrix{Float64}(undef, n_mu, n_T) for col in observable_cols)

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
            mean_N[i, j] = n_avg
            var_N[i, j] = sum(ws .* Ns .^ 2) / sum_w - n_avg^2
            u_avg = sum(ws .* Es) / sum_w
            mean_U[i, j] = u_avg
            var_U[i, j] = sum(ws .* Es .^ 2) / sum_w - u_avg^2
            cov_UN[i, j] = sum(ws .* Es .* Ns) / sum_w - u_avg * n_avg
            for col in observable_cols
                obs_out[col][i, j] = sum(ws .* As[col]) / sum_w
            end
            N_eff[i, j] = sum_w^2 / sum(abs2, ws)
        end
    end

    return (logXi=logXi, mean_N=mean_N, var_N=var_N, mean_U=mean_U, N_eff=N_eff,
            var_U=var_U, cov_UN=cov_UN, observables=obs_out)
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

Return `log(n!)` for non-negative `n`, computed by summing `log(k)`. Accurate
to roundoff accumulation (absolute error ~1e-11 at `n ~ 10^4`), which covers
both the small-`n` atomistic use and the large-`M` `_log_binomial` use.
"""
function _log_factorial(n::Integer)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    n == 0 && return 0.0
    return sum(log, 1:n)
end

"""
    _log_binomial(M::Integer, N::Integer) -> Float64

Return `log(binomial(M, N))` in floating point, computed from `_log_factorial`
so that it does not overflow for lattice sizes where `binomial(M, N)` exceeds
`typemax(Int64)`.
"""
function _log_binomial(M::Integer, N::Integer)
    0 <= N <= M || throw(ArgumentError("require 0 ≤ N ≤ M, got N=$N, M=$M"))
    return _log_factorial(M) - _log_factorial(N) - _log_factorial(M - N)
end

"""
    _fixed_N_log_evidence(df, T_grid; n_walkers, n_cull, ω0, live_energies, kb)
        -> (log_Z_NS::Vector{Float64}, mean_E::Vector{Float64})

Canonical NS evidence `log Z_NS(β) = log Σᵢ ωᵢ exp(−βEᵢ)` and canonical mean
energy at each temperature in `T_grid`, evaluated by log-sum-exp. The discarded
weights `ωᵢ` are built directly in log space, so deep ladders (`iter` beyond
~`700 · n_walkers`, where the linear-space `ωᵢ` underflows Float64) keep their
low-energy samples. When `live_energies` is supplied, the residual prior mass
`ω0 · (K/(K+n_cull))^{n_iters}` (with `n_iters = maximum(df.iter)`) is split
evenly among the supplied energies — the live-set tail correction. An empty
ladder records no compression, so its tail carries the sector's entire prior
mass, exactly 1 and independent of `ω0` (matching the internally-handled
`N = 0` sector). Internal helper shared by the `gc_thermodynamic_stats_fixed_N`
methods.
"""
function _fixed_N_log_evidence(
    df::DataFrame,
    T_grid::AbstractVector{<:Unitful.Temperature};
    n_walkers::Int,
    n_cull::Int,
    ω0::Float64,
    live_energies::Union{Nothing,AbstractVector{<:Real}},
    kb::Float64,
)
    Es = collect(Float64, df.emax)
    if isempty(Es) && (live_energies === nothing || isempty(live_energies))
        throw(ArgumentError(
            "a sector has no discarded samples and no live-set tail; for " *
            "single-configuration sectors supply an empty DataFrame plus a " *
            "live_emax entry carrying the sector energy (e.g. n_walkers " *
            "copies of it)"))
    end
    # Log-space weights: the linear-space ωᵢ underflows Float64 for
    # iter ≳ 700·K, silently dropping the low-energy samples that dominate
    # at low T.
    log_r = log(n_walkers / (n_walkers + n_cull))
    log_ωi = (log(ω0) + log(n_cull / (n_walkers + n_cull))) .+ df.iter .* log_r
    if live_energies === nothing
        Es_all = Es
        log_ws = log_ωi
    else
        Es_live = collect(Float64, live_energies)
        # The live set carries the residual prior mass ω0 · r^n_iters, split
        # evenly among the supplied energies. n_iters comes from the recorded
        # iter values (not the row count), so continuation ladders whose iter
        # does not start at 1 stay consistent with their ωᵢ. An empty ladder
        # records no compression: its tail is the sector's entire prior mass,
        # exactly 1, independent of ω0.
        n_iters = isempty(df.iter) ? 0 : maximum(df.iter)
        log_tail_total = n_iters == 0 ? 0.0 : log(ω0) + n_iters * log_r
        log_tail = log_tail_total - log(length(Es_live))
        Es_all = vcat(Es, Es_live)
        log_ws = vcat(log_ωi, fill(log_tail, length(Es_live)))
    end

    log_Z_NS = Vector{Float64}(undef, length(T_grid))
    mean_E = Vector{Float64}(undef, length(T_grid))
    for (j, T) in enumerate(T_grid)
        β = 1.0 / (kb * ustrip(u"K", T))
        log_terms = log_ws .- β .* Es_all
        max_log = maximum(log_terms)
        ws = exp.(log_terms .- max_log)
        sum_w = sum(ws)
        log_Z_NS[j] = max_log + log(sum_w)
        mean_E[j] = sum(ws .* Es_all) / sum_w
    end
    return log_Z_NS, mean_E
end

"""
    gc_thermodynamic_stats_fixed_N(ns_outputs, N_values, V, atomic_mass, μ_grid, T_grid;
                                   n_walkers=120, n_cull=1, ω0=1.0,
                                   live_emax=nothing, kb=8.617333262e-5)

Atomistic method: compute grand-canonical thermodynamic averages from a stack
of canonical nested-sampling outputs, one per fixed particle number `N`, using
the simulation-box volume `V` and the thermal wavelength. Returns `Xi` in
linear space (contrast the lattice-gas method of this function, which takes
`n_sites` and returns `logXi`).

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
sits in the `K` surviving live walkers. Supplying `live_emax` (one vector of live
walker energies per `N`, normally the K live energies) adds the live-set tail to
each per-N evidence: the residual mass `ω0 · (K/(K+n_cull))^{n_iters}` is split
evenly among the supplied energies. A sector supplied as an empty DataFrame plus
live energies carries mass exactly `1`, independent of `ω0` (see the lattice-gas
method's "Special sectors"). When omitted, the live-set tail is neglected — for
ratio observables (`⟨N⟩`, `⟨U⟩`) the resulting bias is small but visible at low T
or shallow NS; for the absolute `Ξ` it appears as a uniform-in-N prefactor that
does not cancel.

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
    if !allunique(N_values)
        throw(ArgumentError("N_values must not contain duplicate entries"))
    end
    if live_emax !== nothing && length(live_emax) != length(N_values)
        throw(DimensionMismatch("live_emax and N_values must have the same length"))
    end
    if !all(T -> ustrip(u"K", T) > 0, T_grid)
        throw(ArgumentError("every temperature in T_grid must be positive"))
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

        live = live_emax === nothing ? nothing : live_emax[i]
        log_Z_NS[i, :], mean_E_N[i, :] = _fixed_N_log_evidence(
            df, T_grid;
            n_walkers=n_walkers, n_cull=n_cull, ω0=ω0,
            live_energies=live, kb=kb)
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

"""
    gc_thermodynamic_stats_fixed_N(ns_outputs, N_values, n_sites, μ_grid, T_grid;
                                   n_walkers=120, n_cull=1, ω0=1.0,
                                   live_emax=nothing, kb=8.617333262e-5)

Lattice-gas method: compute grand-canonical thermodynamic averages from a
stack of canonical (fixed-`N`) lattice nested-sampling outputs, one per
particle count `N`, on a lattice with `n_sites` sites.

Fixed-`N` lattice NS samples the uniform prior over the `binomial(n_sites, N)`
occupation patterns, so the absolute canonical partition function is
`Z_N(β) = binomial(M, N) · Z_NS^{(N)}(β)` with `M = n_sites` and

```math
\\Xi(\\mu, T) = \\sum_N z^N \\binom{M}{N} Z_{\\mathrm{NS}}^{(N)}(\\beta),
\\qquad z = \\exp(\\beta\\mu).
```

A lattice gas has no momentum integral, so `z = exp(βμ)` directly — no volume
argument and no thermal wavelength enter (contrast the atomistic method of
this function). The sum runs over the supplied `N_values`, which must include
`0`; truncating `N_values` below `n_sites` truncates Ξ accordingly, which is
safe only when `⟨N⟩` at the requested `(μ, T)` sits well below the largest
supplied `N`. (For an athermal hard-core model the truncation at the
close-packing `N` is exact — every higher sector has no allowed
configuration; see the hard-core recipe in the `GenericLatticeHamiltonian`
docstring for the per-`N` sampling setup and the evaluation-temperature
criterion.)

## Special sectors

- `N = 0`: the empty lattice is a single configuration with `E = 0`, so
  `Z_NS^{(0)} = 1`; the corresponding DataFrame contents are ignored.
- `N = n_sites` (and any other single-configuration sector): NS cannot make
  progress because every walker holds the same configuration. Supply an empty
  DataFrame (`DataFrame(iter=Int[], emax=Float64[])`) together with a
  `live_emax` entry holding the sector energy (conventionally `n_walkers`
  copies of it; any count works) — the live-set tail then carries the
  sector's entire prior mass, exactly `1` and independent of `ω0`, matching
  the internally-handled `N = 0` sector. Because such sectors require a
  `live_emax` entry, covering the full range `0:n_sites` in `N_values` is
  only possible with `live_emax` supplied.

## Live-set tail and normalization

As for the atomistic method, the recorded weights `ωᵢ` carry only part of the
prior volume after finite NS termination. Supplying `live_emax` (one vector of
live-walker energies per `N`, normally the `n_walkers` live energies) splits
the residual prior mass `ω0 · (K/(K+n_cull))^{n_iters}` evenly among the
supplied entries of each sector. For fully-normalized absolute `Ξ`, pass
`ω0 = (n_walkers + n_cull)/n_walkers` (Skilling weights) together with
`live_emax`; with the defaults (`ω0 = 1.0`, no tail), `logXi` is biased low,
and the omitted tail additionally leaves an `N`-dependent bias — small but
visible at low `T` or for shallow NS runs, and not exactly cancelling in the
ratio observables `mean_N`, `var_N`, and `mean_U` when per-`N` ladders differ
in length or convergence.

# Arguments
- `ns_outputs::AbstractVector{<:DataFrame}`: one canonical lattice-NS output
  per `N`, each with columns `[:iter, :emax]` (energies in eV, as produced by
  `nested_sampling` on a `LatticeGasWalkers` liveset).
- `N_values::AbstractVector{<:Integer}`: particle counts corresponding to each
  DataFrame. Must include `0`, and every entry must lie in `0:n_sites`.
- `n_sites::Integer`: number of lattice sites `M`.
- `μ_grid::AbstractVector{<:typeof(1.0u"eV")}`: chemical potentials.
- `T_grid::AbstractVector{<:Unitful.Temperature}`: temperatures.

# Keyword arguments
- `n_walkers::Int=120`: number of NS walkers used to produce each DataFrame.
  Must be uniform across `ns_outputs`.
- `n_cull::Int=1`: NS culls per iteration.
- `ω0::Float64=1.0`: initial prior weight for the discarded-sample ladder.
- `live_emax::Union{Nothing,AbstractVector{<:AbstractVector{<:Real}}}=nothing`:
  when supplied, one vector of live-walker energies (in eV) per `N` —
  normally the `n_walkers` live energies; the sector's residual prior mass is
  split evenly among the supplied entries. The entry for `N=0` is ignored.
- `kb::Float64`: Boltzmann constant in eV/K.

# Returns
A `NamedTuple` `(logXi, mean_N, var_N, mean_U, log_Z_N, N_values)`. The first
four fields are `Matrix{Float64}` of size `(length(μ_grid), length(T_grid))`
indexed `[i_μ, i_T]`. `logXi` is the natural log of the absolute grand
partition function — returned in log space (unlike the atomistic method's
`Xi`) because the binomial prior mass grows like `2^M` and `Ξ` overflows
`Float64` for modest lattices. `mean_N` is `⟨N⟩`, `var_N` is `⟨N²⟩ − ⟨N⟩²`,
and `mean_U` is `⟨E⟩` (grand-canonical, in eV). `log_Z_N` is the
`(length(N_values), length(T_grid))` matrix of log canonical
partition-function slices `log[C(M,N) · Z_NS^{(N)}(β)]` — the per-`N`
evidence the assembly is built from — and `N_values` echoes the particle
counts (as `Vector{Int}`) indexing its rows. For an athermal (hard-core)
model evaluated at a sufficiently low `T` (see the hard-core recipe in the
`GenericLatticeHamiltonian` docstring), `log_Z_N` is `log g_N`, the log count
of allowed `N`-particle configurations.
"""
function gc_thermodynamic_stats_fixed_N(
    ns_outputs::AbstractVector{<:DataFrame},
    N_values::AbstractVector{<:Integer},
    n_sites::Integer,
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
        throw(ArgumentError("N_values must include 0 (the empty lattice)"))
    end
    if !all(0 .<= N_values .<= n_sites)
        throw(ArgumentError("every entry of N_values must lie in 0:n_sites"))
    end
    if !allunique(N_values)
        throw(ArgumentError("N_values must not contain duplicate entries"))
    end
    if live_emax !== nothing && length(live_emax) != length(N_values)
        throw(DimensionMismatch("live_emax and N_values must have the same length"))
    end
    if !all(T -> ustrip(u"K", T) > 0, T_grid)
        throw(ArgumentError("every temperature in T_grid must be positive"))
    end

    N_int = collect(Int, N_values)
    n_N = length(N_int)
    n_mu = length(μ_grid)
    n_T = length(T_grid)

    log_Z_NS = Matrix{Float64}(undef, n_N, n_T)
    mean_E_N = Matrix{Float64}(undef, n_N, n_T)

    for (i, df) in enumerate(ns_outputs)
        # The N = 0 sector is the empty lattice: a single configuration with
        # E = 0, so Z_NS^{(0)} = 1. The corresponding DataFrame is ignored.
        if N_int[i] == 0
            log_Z_NS[i, :] .= 0.0
            mean_E_N[i, :] .= 0.0
            continue
        end

        live = live_emax === nothing ? nothing : live_emax[i]
        log_Z_NS[i, :], mean_E_N[i, :] = _fixed_N_log_evidence(
            df, T_grid;
            n_walkers=n_walkers, n_cull=n_cull, ω0=ω0,
            live_energies=live, kb=kb)
    end

    log_binom = [_log_binomial(n_sites, N) for N in N_int]

    logXi = Matrix{Float64}(undef, n_mu, n_T)
    mean_N = Matrix{Float64}(undef, n_mu, n_T)
    var_N = Matrix{Float64}(undef, n_mu, n_T)
    mean_U = Matrix{Float64}(undef, n_mu, n_T)

    for (j, T) in enumerate(T_grid)
        β = 1.0 / (kb * ustrip(u"K", T))
        for (k, μ) in enumerate(μ_grid)
            μ_val = ustrip(u"eV", μ)

            log_w = log_Z_NS[:, j] .+ N_int .* (β * μ_val) .+ log_binom
            max_log = maximum(log_w)
            ws = exp.(log_w .- max_log)
            sum_w = sum(ws)

            logXi[k, j] = max_log + log(sum_w)
            mean_N[k, j] = sum(ws .* N_int) / sum_w
            mean_N2 = sum(ws .* (N_int .^ 2)) / sum_w
            var_N[k, j] = mean_N2 - mean_N[k, j]^2
            mean_U[k, j] = sum(ws .* view(mean_E_N, :, j)) / sum_w
        end
    end

    return (logXi=logXi, mean_N=mean_N, var_N=var_N, mean_U=mean_U,
            log_Z_N=log_Z_NS .+ log_binom, N_values=N_int)
end


include("microcanonical_inflection.jl")

end # module AnalysisTools