# =============================================================================
# Microcanonical inflection-point analysis of nested-sampling output.
#
# Implements the microcanonical entropy and its energetic derivatives following
# Schnabel et al., Phys. Rev. E 84, 011127 (2011) and Qi & Bachmann, Phys. Rev.
# Lett. 120, 180601 (2018). Transitions appear as inflection points of the
# inverse microcanonical temperature β(E)=dS/dE; the order is read from the
# higher derivatives (see `inflection_transitions`).
#
# Nested sampling gives the entropy essentially exactly along the (contiguous)
# cull index i: the enclosed prior volume is X_i = (K/(K+1))^i, so
# ln X_i = i·ln(K/(K+1)). With the energy ladder E(i) (E = `emax`), the caloric
# entropy is S(E) = ln g(E) = i·ln(K/(K+1)) − ln|dE/di| (+const), and derivatives
# are taken against the dense, uniform i-axis (local cubic fits) — far more stable
# than differentiating a binned ln g in energy space.
#
# NUMERICS NOTE: γ=d²S/dE² amplifies the finite-walker "staircase" granularity of
# the NS ladder (step ~ 1/(K·g)). That roughness shrinks with the NUMBER OF
# WALKERS K (not with more independent runs — the ladder is reproducible across
# seeds). Smoothing here is a light, necessary differentiation aid, NOT a
# substitute for adequate K: heavy smoothing on under-converged ladders biases the
# transition temperatures. Use `transition_convergence` to check K-convergence.
# =============================================================================

# Local cubic least-squares fit of `y` vs `x` over the window [idx-W, idx+W];
# returns the value, first, and second derivative at node `idx`.
function _local_cubic(x::AbstractVector, y::AbstractVector, idx::Int, W::Int)
    lo = max(firstindex(x), idx - W)
    hi = min(lastindex(x), idx + W)
    dx = x[lo:hi] .- x[idx]
    V = [dx[k]^p for k in eachindex(dx), p in 0:3]
    c = V \ y[lo:hi]
    return c[1], c[2], 2 * c[3]
end

# First derivative of `y` w.r.t. `x` at every point, via local cubic fits.
function _derivative(x::AbstractVector, y::AbstractVector; W::Int = 12)
    d = similar(float.(y))
    @inbounds for j in eachindex(x)
        lo = max(firstindex(x), j - W)
        hi = min(lastindex(x), j + W)
        dx = x[lo:hi] .- x[j]
        V = [dx[k]^p for k in eachindex(dx), p in 0:3]
        d[j] = (V \ y[lo:hi])[2]
    end
    return d
end

# Core: build the smoothed microcanonical curves on a set of i-space nodes.
# Returns a NamedTuple with energies `E` (ascending) and, as requested,
# `S_vol`, `S_caloric`, `β`, `γ`, `δ`.
function _microcanonical_core(df::DataFrame, n_walkers::Int; n_cull::Int = 1,
                              n_nodes::Int = 400, halfwidth::Int = 0, max_order::Int = 2)
    n_walkers > 0 || throw(ArgumentError("n_walkers must be positive"))
    1 <= max_order <= 3 || throw(ArgumentError("max_order must be 1, 2, or 3"))
    hasproperty(df, :iter) && hasproperty(df, :emax) ||
        throw(ArgumentError("df must have :iter and :emax columns (nested_sampling output)"))
    i = Float64.(df.iter)
    e = Float64.(df.emax)
    n = length(i)
    n >= 9 || throw(ArgumentError("NS ladder too short (need ≥ 9 recorded iterations)"))

    cc = log(n_walkers / (n_walkers + n_cull))         # ln(K/(K+n_cull)) < 0
    W = halfwidth > 0 ? halfwidth : clamp(n ÷ 55, 25, 2000)
    W = min(W, (n - 3) ÷ 2)
    nodes = unique(round.(Int, range(W + 1, n - W; length = min(n_nodes, n - 2W))))

    E = Float64[]; Svol = Float64[]; Scal = Float64[]; β = Float64[]; slope = Float64[]
    for idx in nodes
        _, ep, epp = _local_cubic(i, e, idx, W)        # dE/di, d²E/di²
        ep < 0 || continue                              # ladder must descend
        push!(E, e[idx]); push!(slope, ep)
        push!(Svol, cc * i[idx])                        # volume entropy ln G = ln X
        push!(Scal, cc * i[idx] - log(-ep))             # caloric entropy ln g (+const)
        push!(β, cc / ep - epp / ep^2)                  # β = dS/dE (caloric)
    end
    length(E) >= 5 || throw(ArgumentError("too few usable nodes; check the ladder/halfwidth"))

    o = sortperm(E)                                     # ascending in energy
    E, Svol, Scal, β, slope = E[o], Svol[o], Scal[o], β[o], slope[o]
    γ = max_order >= 2 ? _derivative(E, β) : Float64[]
    δ = max_order >= 3 ? _derivative(E, γ) : Float64[]
    return (E = E, S_vol = Svol, S_caloric = Scal, β = β, γ = γ, δ = δ, slope = slope)
end

"""
    microcanonical_entropy(df::DataFrame, n_walkers::Int;
                           n_cull=1, kind=:caloric, n_nodes=400, halfwidth=0)

Estimate the microcanonical entropy ``S(E)`` from a canonical nested-sampling
output `df` (columns `:iter`, `:emax`) produced with `n_walkers` live points.

Along the contiguous cull index ``i`` the enclosed prior volume is
``X_i=(K/(K+n_{cull}))^i``, so the volume entropy is ``S_{vol}=\\ln X_i=i\\ln(K/(K+n_{cull}))``
and the caloric entropy is ``S(E)=\\ln g(E)=i\\ln(K/(K+n_{cull}))-\\ln|dE/di|``,
with the ladder slope ``dE/di`` obtained from local cubic fits against the dense,
uniform `i`-axis.

# Arguments
- `kind`: `:caloric` for ``S=\\ln g`` (matches Schnabel/Qi–Bachmann; default) or
  `:volume` for the Hertz/Gibbs volume entropy ``S_{vol}=\\ln G`` (one fewer
  derivative, slightly cleaner).
- `n_cull`: NS culls per iteration (default 1).
- `n_nodes`: number of smoothed output nodes (default 400).
- `halfwidth`: half-window (in iterations) of the local cubic smoother; `0` (default)
  picks `clamp(n÷55, 25, 2000)`.

# Returns
`(E, S)` — energies (ascending, same units as `:emax`) and entropy (dimensionless,
up to an additive constant).

See also [`caloric_derivatives`](@ref), [`inflection_transitions`](@ref).
"""
function microcanonical_entropy(df::DataFrame, n_walkers::Int; n_cull::Int = 1,
                                kind::Symbol = :caloric, n_nodes::Int = 400, halfwidth::Int = 0)
    kind in (:caloric, :volume) || throw(ArgumentError("kind must be :caloric or :volume"))
    r = _microcanonical_core(df, n_walkers; n_cull = n_cull, n_nodes = n_nodes,
                             halfwidth = halfwidth, max_order = 1)
    return r.E, kind === :caloric ? r.S_caloric : r.S_vol
end

"""
    caloric_derivatives(df::DataFrame, n_walkers::Int;
                        n_cull=1, max_order=2, n_nodes=400, halfwidth=0)

Microcanonical entropy derivatives from a nested-sampling output `df`:
``β(E)=dS/dE`` (inverse caloric temperature; the microcanonical temperature is
``k_BT=1/β``), ``γ(E)=d^2S/dE^2``, and optionally ``δ(E)=d^3S/dE^3``.

Derivatives are computed in iteration-index space (see [`microcanonical_entropy`])
and then differentiated against energy with local cubic fits.

# Arguments
- `max_order`: highest derivative to return (1 → β only; 2 → β, γ; 3 → β, γ, δ).
- `n_cull`, `n_nodes`, `halfwidth`: as in [`microcanonical_entropy`].

# Returns
A `NamedTuple` with `E` (ascending energies, same units as `:emax`) and `β` (units
of inverse energy), plus `γ` and `δ` as requested by `max_order`. Note `γ`, `δ` are
sensitive to NS statistics — increase the walker count `K` (not the number of seeds)
for cleaner higher derivatives; see [`transition_convergence`](@ref).
"""
function caloric_derivatives(df::DataFrame, n_walkers::Int; n_cull::Int = 1,
                             max_order::Int = 2, n_nodes::Int = 400, halfwidth::Int = 0)
    1 <= max_order <= 3 || throw(ArgumentError("max_order must be 1, 2, or 3"))
    r = _microcanonical_core(df, n_walkers; n_cull = n_cull, n_nodes = n_nodes,
                             halfwidth = halfwidth, max_order = max_order)
    max_order == 1 && return (E = r.E, β = r.β)
    max_order == 2 && return (E = r.E, β = r.β, γ = r.γ)
    return (E = r.E, β = r.β, γ = r.γ, δ = r.δ)
end

# median / linear-interpolated quantile without a Statistics dependency
function _median(v::AbstractVector)
    s = sort(v); m = length(s)
    m == 0 && return zero(eltype(s))
    iseven(m) ? (s[m÷2] + s[m÷2+1]) / 2 : s[m÷2+1]
end
function _quantile(v::AbstractVector, q::Real)
    s = sort(collect(v)); m = length(s)
    m == 1 && return float(s[1])
    h = (m - 1) * q + 1
    lo = floor(Int, h); hi = min(lo + 1, m)
    return s[lo] + (h - lo) * (s[hi] - s[lo])
end

# Local extrema of the order-n entropy derivative `Sn` carrying the Qi–Bachmann
# transition signature: odd order → positive local minimum, even order → negative
# local maximum. Prominence-filtered (relative to the interior range) and merged by
# a minimum energy separation. Returns node indices into `E`/`Sn`.
function _order_extrema(E::Vector{Float64}, Sn::Vector{Float64}, order::Int;
                        prominence::Float64, edge::Int, min_sep::Float64)
    n = length(Sn)
    lo, hi = 1 + edge, n - edge
    lo < hi || return Int[]
    interior = view(Sn, lo:hi)
    rng = _quantile(interior, 0.9) - _quantile(interior, 0.1)   # robust to divergence outliers
    rng > 0 || return Int[]
    thr = prominence * rng
    odd = isodd(order)
    cands = Tuple{Int,Float64}[]
    for i in lo:hi
        if odd                                            # positive local minimum
            (Sn[i] < Sn[i-1] && Sn[i] <= Sn[i+1] && Sn[i] > 0) || continue
            l = i; while l > 1 && Sn[l-1] >= Sn[i]; l -= 1; end
            r = i; while r < n && Sn[r+1] >= Sn[i]; r += 1; end
            prom = min(maximum(view(Sn, l:i)), maximum(view(Sn, i:r))) - Sn[i]
        else                                              # negative local maximum
            (Sn[i] > Sn[i-1] && Sn[i] >= Sn[i+1] && Sn[i] < 0) || continue
            l = i; while l > 1 && Sn[l-1] <= Sn[i]; l -= 1; end
            r = i; while r < n && Sn[r+1] <= Sn[i]; r += 1; end
            prom = Sn[i] - max(minimum(view(Sn, l:i)), minimum(view(Sn, i:r)))
        end
        prom >= thr && push!(cands, (i, prom))
    end
    sort!(cands, by = c -> -c[2])
    chosen = Int[]
    for (i, _) in cands
        all(abs(E[i] - E[j]) > min_sep for j in chosen) && push!(chosen, i)
    end
    return sort(chosen)
end

"""
    inflection_transitions(df::DataFrame, n_walkers::Int;
        n_cull=1, max_order=2, n_nodes=400, halfwidth=0,
        kb=8.617333262e-5, prominence=0.05, edge=4,
        min_separation=nothing, energy_window=nothing, beta_max=nothing)

Identify and classify phase transitions in a nested-sampling output `df` by
microcanonical inflection-point analysis (Schnabel et al. 2011; Qi & Bachmann
2018). Returns one entry per transition, ordered by energy.

A transition of order `n` is a "least-sensitive" inflection point of the entropy,
detected as an extremum of the `n`-th derivative `Sⁿ(E)` with the Qi–Bachmann sign:
**odd order → a positive local minimum** of `Sⁿ` (order 1 = β backbending,
first-order), **even order → a negative local maximum** of `Sⁿ` (order 2 = γ peak,
second-order). Orders `1…max_order` are searched. The transition temperature is
`T_tr = 1/(kb·β(E_tr))`.

# Arguments
- `max_order`: highest transition order to search (1–3).
- `kb`: Boltzmann constant; default `8.617333262e-5` eV/K gives `T_tr` in Kelvin
  when `:emax` is in eV. Pass `kb=1.0` for the microcanonical temperature in energy
  units (`kT`).
- `prominence`: peak prominence threshold as a fraction of the (robust 10–90
  percentile) derivative range (default 0.20). Raise it to keep only strong
  transitions.
- `ground_trim`: fraction of the energy span (above the ground state) to discard
  before differentiating, removing the `β → ∞` divergence as `dE/di → 0` (default
  0.05). γ/δ are recomputed on the trimmed window so transitions are not
  contaminated by the divergence.
- `edge`, `min_separation`, `energy_window`, `beta_max`: numerical controls;
  `min_separation` defaults to 3% of the (trimmed) energy span, and `beta_max`
  optionally adds an explicit β ceiling.
- `n_cull`, `n_nodes`, `halfwidth`: as in [`caloric_derivatives`].

# Returns
`Vector{NamedTuple}` with fields `E_tr`, `T_tr`, `order::Int`, `kind`
(`:independent` or `:dependent`), and `strength` (the `Sⁿ` extremum value).
`kind` uses the Qi–Bachmann heuristic: a transition accompanied by a lower-order
transition at lower energy is `:dependent`, otherwise `:independent`.

!!! note
    γ and higher derivatives are sensitive to NS statistics. The ladder roughness
    is finite-walker granularity (it does **not** average out over seeds), so use
    enough walkers `K`; verify with [`transition_convergence`](@ref). Heavy
    smoothing (large `halfwidth`) on under-converged ladders biases `T_tr`.
"""
function inflection_transitions(df::DataFrame, n_walkers::Int; n_cull::Int = 1,
        max_order::Int = 2, n_nodes::Int = 400, halfwidth::Int = 0,
        kb::Float64 = 8.617333262e-5, prominence::Float64 = 0.20, edge::Int = 4,
        ground_trim::Float64 = 0.05,
        min_separation::Union{Nothing,Float64} = nothing,
        energy_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
        beta_max::Union{Nothing,Float64} = nothing)
    1 <= max_order <= 3 || throw(ArgumentError("max_order must be 1, 2, or 3"))
    0 <= ground_trim < 0.5 || throw(ArgumentError("ground_trim must be in [0, 0.5)"))
    r = _microcanonical_core(df, n_walkers; n_cull = n_cull, n_nodes = n_nodes,
                             halfwidth = halfwidth, max_order = 1)
    E, β = r.E, r.β
    # Trim the ground-state β divergence (the ladder flattens, dE/di → 0, β → ∞)
    # BEFORE differentiating, so γ/δ near the transitions are not contaminated.
    Emin, Emax = extrema(E)
    Elo = Emin + ground_trim * (Emax - Emin)
    keep = (β .> 0) .& (E .>= Elo)
    beta_max === nothing || (keep = keep .& (β .< beta_max))
    if energy_window !== nothing
        keep = keep .& (E .>= energy_window[1]) .& (E .<= energy_window[2])
    end
    idx = findall(keep)
    length(idx) >= 2edge + 5 || return NamedTuple[]
    Ew, βw = E[idx], β[idx]
    # derivatives recomputed on the windowed β (free of the ground divergence)
    γw = max_order >= 2 ? _derivative(Ew, βw) : Float64[]
    δw = max_order >= 3 ? _derivative(Ew, γw) : Float64[]
    derivs = (βw, γw, δw)
    msep = min_separation === nothing ? 0.03 * (maximum(Ew) - minimum(Ew)) : min_separation
    found = NamedTuple[]
    for n in 1:max_order
        Sn = derivs[n]
        for i in _order_extrema(Ew, Sn, n; prominence = prominence, edge = edge, min_sep = msep)
            push!(found, (E_tr = Ew[i], T_tr = 1 / (kb * βw[i]), order = n, strength = Sn[i]))
        end
    end
    sort!(found, by = t -> t.E_tr)
    return [(E_tr = t.E_tr, T_tr = t.T_tr, order = t.order,
             kind = any(s.order < t.order && s.E_tr < t.E_tr for s in found) ? :dependent : :independent,
             strength = t.strength) for t in found]
end

# closest same-order transition to `t` within `match_tol` (relative in T_tr); else nothing
function _closest(transitions, t, match_tol::Float64)
    best = nothing; bestd = Inf
    for s in transitions
        s.order == t.order || continue
        d = abs(s.T_tr - t.T_tr)
        if d <= match_tol * t.T_tr && d < bestd
            best, bestd = s, d
        end
    end
    return best
end

"""
    transition_convergence(dfs, n_walkers; tol=0.1, match_tol=0.25, kwargs...)

Assess walker-count (`K`) convergence of microcanonical inflection-point
transitions. `dfs` is a collection of nested-sampling outputs at the walker counts
`n_walkers` (any order; they are sorted ascending internally). `inflection_transitions`
is run at each `K` and the results from the **largest `K`** are taken as the current
best estimate, then traced down the `K`-ladder.

This is the recommended way to use the inflection analysis: the γ (and higher)
derivatives carry finite-`K` "staircase" granularity that does **not** average out
over independent runs, so a transition should only be trusted once its temperature
and order stop drifting as `K` increases. Near-ground transitions converge last.

# Arguments
- `tol`: relative `T_tr` drift between the two largest `K` below which a transition
  is flagged `converged` (default 0.1).
- `match_tol`: relative `T_tr` window for matching the same transition across `K`
  (default 0.25).
- `kwargs...`: forwarded to [`inflection_transitions`] (e.g. `max_order`, `prominence`,
  `ground_trim`, `kb`).

# Returns
`Vector{NamedTuple}`, one per transition found at the largest `K`, with fields
`T_tr`, `order`, `kind`, `converged::Bool`, `T_drift` (relative change from the
second-largest `K`, `Inf` if it has no match there), `T_by_K` (the matched `T_tr`
at each `K`, `missing` where unmatched), and `n_walkers` (the sorted `K` values).
"""
function transition_convergence(dfs, n_walkers; tol::Float64 = 0.1,
                                match_tol::Float64 = 0.25, kwargs...)
    length(dfs) == length(n_walkers) ||
        throw(DimensionMismatch("dfs and n_walkers must have the same length"))
    length(dfs) >= 2 ||
        throw(ArgumentError("need at least two walker counts to assess convergence"))
    perm = sortperm(collect(n_walkers))
    dfv = collect(dfs)[perm]
    Ks = collect(n_walkers)[perm]
    per_K = [inflection_transitions(dfv[k], Ks[k]; kwargs...) for k in eachindex(dfv)]
    ref = per_K[end]
    out = NamedTuple[]
    for t in ref
        T_by_K = Union{Missing,Float64}[
            (c = _closest(per_K[k], t, match_tol); c === nothing ? missing : c.T_tr)
            for k in eachindex(per_K)]
        prevc = _closest(per_K[end-1], t, match_tol)
        drift = prevc === nothing ? Inf : abs(prevc.T_tr - t.T_tr) / t.T_tr
        push!(out, (T_tr = t.T_tr, order = t.order, kind = t.kind,
                    converged = drift <= tol, T_drift = drift,
                    T_by_K = T_by_K, n_walkers = Ks))
    end
    return out
end
