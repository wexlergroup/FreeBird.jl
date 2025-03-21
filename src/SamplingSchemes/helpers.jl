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