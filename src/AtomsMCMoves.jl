module AtomsMCMoves

using ExtXYZ
# using LinearAlgebra
using Setfield
using Distributions

using ..Potentials
using ..AbstractWalkers

export random_walk


"""
    random_walk(n_steps::Int, at::Atoms, lj::LJParameters, step_size::Float64, emax::Float64, frozen::Int)
Perform a random walk of `n_steps` steps on the atoms in `at` using a step size of `step_size`.
"""
function random_walk(n_steps::Int, at::Atoms, lj::LJParameters, step_size::Float64, emax::Float64, frozen::Int)
    n_accept = 0
    for i_mc_step in 1:n_steps
        for i_at in (frozen+1):length(at)
            dx = rand(Uniform(-step_size,step_size))
            dy = rand(Uniform(-step_size,step_size))
            dz = rand(Uniform(-step_size,step_size))
            orig_energy = at.system_data.energy
            # println("orig_energy: ", orig_energy) # debug
            orig_pos = at.atom_data.position[i_at]
            # println("orig_pos: ", orig_pos) # debug
            at = @set at.atom_data.position[i_at][1].val += dx
            at = @set at.atom_data.position[i_at][2].val += dy
            at = @set at.atom_data.position[i_at][3].val += dz
            energy = compute_total_energy(at, lj)
            if energy >= emax
                # println("energy too high: ", energy, " >= ", emax) # debug
                at = @set at.atom_data.position[i_at] = orig_pos
            else
                at = @set at.system_data.energy = energy
                # println("accepted energy: ", energy) # debug
                n_accept += 1
            end
        end
    end
    # println("after walk at.system_data.energy: ", at.system_data.energy) # debug
    return n_accept, at
end



end # module AtomsMonteCarloMoves