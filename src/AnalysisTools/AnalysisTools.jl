"""
    AnalysisTools

Module for analyzing the output of the sampling.
"""
module AnalysisTools

using DataFrames
using CSV, Arrow

export read_output
export ωᵢ, partition_function, internal_energy, cv
export gc_thermodynamic_stats
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

    log_prior_mass = n_sites * log1p(z0)
    log_z0 = log(z0)

    n_mu = length(μs)
    n_T = length(Ts)
    logXi = Matrix{Float64}(undef, n_mu, n_T)
    mean_N = Matrix{Float64}(undef, n_mu, n_T)
    var_N = Matrix{Float64}(undef, n_mu, n_T)
    mean_U = Matrix{Float64}(undef, n_mu, n_T)
    N_eff = Matrix{Float64}(undef, n_mu, n_T)

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
            mean_U[i, j] = sum(ws .* Es) / sum_w
            N_eff[i, j] = sum_w^2 / sum(abs2, ws)
        end
    end

    return (logXi=logXi, mean_N=mean_N, var_N=var_N, mean_U=mean_U, N_eff=N_eff)
end


include("microcanonical_inflection.jl")

end # module AnalysisTools