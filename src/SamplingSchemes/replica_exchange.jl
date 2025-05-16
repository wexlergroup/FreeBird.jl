"""
    ReplicaExchangeParameters

A structure to hold the parameters for the replica exchange sampling scheme.

# Fields
- `temperatures::Vector{Float64}`: The temperatures of the replicas in ascending order.
- `equilibrium_steps::Int64`: The number of steps to equilibrate each replica.
- `sampling_steps::Int64`: The number of steps to sample each replica.
- `swap_interval::Int64`: The number of steps between attempted replica-exchange moves.
- `random_seed::Int64`: The seed for the random number generator.
"""
mutable struct ReplicaExchangeParameters <: SamplingParameters
    temperatures::Vector{Float64}
    equilibrium_steps::Int64
    sampling_steps::Int64
    swap_interval::Int64
    random_seed::Int64
    function ReplicaExchangeParameters(
        temperatures;
        equilibrium_steps::Int64=10_000,
        sampling_steps::Int64=10_000,
        swap_interval::Int64=100,
        random_seed::Int64=1234
    )
        new(temperatures, equilibrium_steps, sampling_steps, swap_interval, random_seed)
    end
end
