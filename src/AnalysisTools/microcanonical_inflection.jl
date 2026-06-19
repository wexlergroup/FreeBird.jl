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

    E = Float64[]; Svol = Float64[]; Scal = Float64[]; β = Float64[]
    for idx in nodes
        _, ep, epp = _local_cubic(i, e, idx, W)        # dE/di, d²E/di²
        ep < 0 || continue                              # ladder must descend
        push!(E, e[idx])
        push!(Svol, cc * i[idx])                        # volume entropy ln G = ln X
        push!(Scal, cc * i[idx] - log(-ep))             # caloric entropy ln g (+const)
        push!(β, cc / ep - epp / ep^2)                  # β = dS/dE (caloric)
    end
    length(E) >= 5 || throw(ArgumentError("too few usable nodes; check the ladder/halfwidth"))

    o = sortperm(E)                                     # ascending in energy
    E, Svol, Scal, β = E[o], Svol[o], Scal[o], β[o]
    γ = max_order >= 2 ? _derivative(E, β) : Float64[]
    δ = max_order >= 3 ? _derivative(E, γ) : Float64[]
    return (E = E, S_vol = Svol, S_caloric = Scal, β = β, γ = γ, δ = δ)
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
