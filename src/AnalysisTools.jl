module AnalysisTools

using DataFrames
using CSV

export read_df
export gamma_factors, partition_function, internal_energy, cv

function read_df(filename)
    return DataFrame(CSV.File(filename))
end

function gamma_factors(iters::Vector{Int}, n_walkers::Int)
    gi = (1/(n_walkers+1))*(n_walker/(n_walkers+1)).^iters
    return gi
end

function partition_function(beta::Float64, 
                                    gi::Vector{Float64}, 
                                    ei::Vector{Float64})
    z = sum(gi.*exp.(-ei.*beta))
    return z
end

function internal_energy(beta::Float64, 
                                 gi::Vector{Float64}, 
                                 ei::Vector{Float64})
    u = sum(gi.*ei.*exp.(-ei.*beta))/sum(gi.*exp.(-ei.*beta))
    return u
end

function cv(beta::Float64, 
                                       gi::Vector{Float64}, 
                                       ei::Vector{Float64}
                                       dof::Int,)
    z = partition_function(beta, gi, ei)
    u = internal_energy(beta, gi, ei)
    kb = 8.617333262e-5 # eV/K
    cv = dof*kb/2.0 + kb*beta^2 * (sum(gi.*ei.^2 .*exp.(-ei.*beta))/z - u^2)
    return cv
end


end # module AnalysisTools