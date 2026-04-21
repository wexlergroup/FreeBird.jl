"""
    adjust_step_size(params::SamplingParameters, rate::Float64)

Adjusts the step size of the sampling algorithm based on the acceptance rate. 
The step size is increased by 10% if the acceptance rate is above the upper limit of the range,
and decreased by 10% if the acceptance rate is below the lower limit of the range.

# Arguments
- `params::SamplingParameters`: The parameters of the sampling algorithm.
- `rate::Float64`: The acceptance rate of the algorithm.
- `range::Tuple{Float64, Float64}`: The range of acceptance rates for adjusting the step size. Default is (0.25, 0.75).

# Returns
- `params::SamplingParameters`: The updated parameters with adjusted step size.
"""
function adjust_step_size(params::SamplingParameters, rate::Float64; range::Tuple{Float64, Float64}=(0.25, 0.75))
    if rate > range[2] && params.step_size < params.step_size_up
        params.step_size *= 1.1
    elseif rate < range[1] && params.step_size > params.step_size_lo
        params.step_size *= 0.9
    end
    return params
end

"""
    adjust_cluster_p(params::SamplingParameters, rate::Float64, iteration::Int; target::Float64=0.3, floor::Float64=0.01, ceiling::Float64=1.0)

Adjusts the cluster growth probability based on the acceptance rate.
Uses a simple multiplicative rule: `p *= 0.9` if rate is below target,
`p *= 1.1` if rate is at or above target, clamped to `[floor, ceiling]`.

Also logs the adjusted `cluster_p`, acceptance rate, and NS iteration index
to `params.cluster_p_history`, `params.cluster_accept_history`, and
`params.cluster_adjust_iterations` for post-run diagnostics.

# Arguments
- `params::SamplingParameters`: The parameters containing `cluster_p`.
- `rate::Float64`: The cluster move acceptance rate for the current window.
- `iteration::Int`: The current NS iteration index.
- `target::Float64`: The target acceptance rate (default 0.3).
- `floor::Float64`: Lower bound for cluster_p (default 0.01).
- `ceiling::Float64`: Upper bound for cluster_p (default 1.0).

# Returns
- `params::SamplingParameters`: The updated parameters with adjusted cluster_p.
"""
function adjust_cluster_p(params::SamplingParameters, rate::Float64, iteration::Int; target::Float64=0.3, floor::Float64=0.01, ceiling::Float64=1.0)
    if rate < target
        params.cluster_p *= 0.9
    else
        params.cluster_p *= 1.1
    end
    params.cluster_p = clamp(params.cluster_p, floor, ceiling)
    push!(params.cluster_p_history, params.cluster_p)
    push!(params.cluster_accept_history, rate)
    push!(params.cluster_adjust_iterations, iteration)
    return params
end