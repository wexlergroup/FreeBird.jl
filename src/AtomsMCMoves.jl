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
function random_walk(n_steps::Int, at::Atoms, lj::LJParameters, step_size::Float64, emax::Float64; 
                    frozen::Int=0, e_shift::Float64=0.0)
    n_accept = 0
    for i_mc_step in 1:n_steps
        for i_at in (frozen+1):length(at)
            dx, dy, dz = [rand(Uniform(-step_size,step_size)) for _ in 1:3]
            # println("orig_energy: ", at.system_data.energy) # debug
            orig_pos = deepcopy(at.atom_data.position[i_at])
            # println("orig_pos: ", orig_pos) # debug
            at = @set at.atom_data.position[i_at][1].val += dx
            at = @set at.atom_data.position[i_at][2].val += dy
            at = @set at.atom_data.position[i_at][3].val += dz
            energy = interaction_energy(at, lj; frozen=frozen) + e_shift
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




function single_atom_demon_walk(at::Atoms, lj::LJParameters, step_size::Float64; 
                    frozen::Int=0, e_shift::Float64=0.0, e_demon::Float64=0.0)
    accept = false
    i_at = rand((frozen+1):length(at))
    orig_pos = deepcopy(at.atom_data.position[i_at])
    orig_energy = at.system_data.energy
    dx, dy, dz = [rand(Uniform(-step_size,step_size)) for _ in 1:3]
    at = @set at.atom_data.position[i_at][1].val += dx
    at = @set at.atom_data.position[i_at][2].val += dy
    at = @set at.atom_data.position[i_at][3].val += dz
    new_energy = interaction_energy(at, lj; frozen=frozen) + e_shift
    ΔE = new_energy - orig_energy
    if ΔE <= 0
        e_demon += ΔE
        accept = true
        at = @set at.system_data.energy = energy
    elseif 0 < ΔE <= e_demon
        e_demon -= ΔE
        accept = true
        at = @set at.system_data.energy = energy
    else
        at = @set at.atom_data.position[i_at] = orig_pos
    end
    return accept, at, e_demon
end

function MC_nve_walk(n_steps::Int, at::Atoms, lj::LJParameters, step_size::Float64, emax::Float64; 
    frozen::Int=0, e_shift::Float64=0.0)
    e_demon = 0.0
    accempt_count = 0
    for _ in 1:n_steps
        accept, at, e_demon = single_atom_demon_walk(at, lj, step_size; frozen=frozen, e_shift=e_shift, e_demon=e_demon)
        if accept
            accempt_count += 1
        end
    end
    return accempt_count/n_steps, at
end


end # module AtomsMonteCarloMoves