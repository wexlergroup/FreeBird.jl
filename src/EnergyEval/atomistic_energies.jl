"""
    pbc_dist(pos1, pos2, at)

Compute the distance between two positions considering periodic boundary conditions. Currently only works for orthorhombic lattices.

# Arguments
- `pos1::Union{SVector{T},Vector{T}}`: The first position.
- `pos2::Union{SVector{T},Vector{T}}`: The second position.
- `at::AbstractSystem`: The abstract system containing boundary conditions and bounding box.

# Returns
- `d::Float64`: The distance between `pos1` and `pos2` considering periodic boundary conditions.

"""
function pbc_dist(pos1::Union{SVector{T},Vector{T}},
                  pos2::Union{SVector{T},Vector{T}},
                  at::AbstractSystem) where {T}
    pbc = periodicity(at)
    box = cell_vectors(at)
    distsq = 0.0u"Å"^2
    for i in eachindex(pos1)
        if pbc[i] == true
            distsq += min(abs(pos1[i] - pos2[i]), box[i][i] - abs(pos1[i] - pos2[i]))^2
        elseif pbc[i] == false
            distsq += (pos1[i] - pos2[i])^2
        end
    end
    return sqrt(distsq)
end



"""
    pair_energy(r::typeof(1.0u"Å"), lj::LJParameters)
Compute the energy of a pair of particles separated by distance `r` using the Lennard-Jones potential.
# Arguments
- `r::typeof(1.0u"Å")`: The distance between the two particles.
- `lj::LJParameters`: The Lennard-Jones parameters.
# Returns
- `energy::typeof(0.0u"eV")`: The energy of the pair of particles.

"""
pair_energy(r::typeof(1.0u"Å"), lj::LJParameters) = lj_energy(r, lj)


"""
    inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::AbstractPotential)

Compute the energy between two components of a system using a specified (pairwise) potential.

# Arguments
- `at1::AbstractSystem`: The first component of the system.
- `at2::AbstractSystem`: The second component of the system.
- `pot::AbstractPotential`: The potential used to compute the energy.

# Returns
- `energy`: The energy between the two components.

"""
function inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::AbstractPotential)
    # build pairs of particles
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2)]
    # @show pairs # DEBUG
    energy = Array{typeof(0.0u"eV"), 1}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        # @show i,j # DEBUG
        (i, j) = pairs[k]
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energy[k] = pair_energy(r, pot)
    end
    # energy = energy*u"eV"
    return sum(energy)
end

"""
    intra_component_energy(at::AbstractSystem, pot::AbstractPotential)

Compute the energy within a component of a system using a specified (pairwise) potential.

# Arguments
- `at::AbstractSystem`: The component of the system.
- `pot::AbstractPotential`: The potential used to compute the energy.

# Returns
- `energy`: The energy within the component.

"""

#testing stuff that is added on -- delete later
function debug_imports()
    println("Testing imports in atomistic_energies.jl:")
    println("SMatrix: ", @isdefined SMatrix)
    println("unitcell: ", @isdefined unitcell)
    println("ustrip: ", @isdefined ustrip)
    println("position: ", @isdefined position)
end

# === GPU-Compatible Struct for Lennard-Jones Parameters ===
struct LJParams
    epsilon::Float64
    sigma::Float64
end

# === GPU Kernel: Compute Pair Energies (Intra-component) ===
function compute_intra_pair_energies_kernel!(
    energy::CuDeviceVector{Float64},
    pairs::CuDeviceVector{NTuple{2, Int}},
    pos::CuDeviceVector{SVector{3, Float64}},
    pot_params::LJParams,
    cell::SMatrix{3,3,Float64}
)
    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if k >= 1 && k <= length(pairs)
        (i, j) = pairs[k]
        r_vec = pbc_dist_gpu(pos[i], pos[j], cell)
        energy[k] = pair_energy_gpu(norm(r_vec), pot_params)
    end
    return
end

# === Host Function: GPU version of `intra_component_energy` ===
function intra_component_energy(at::AbstractSystem, pot::AbstractPotential)
    # Early return for small systems
    if length(at) < 2
        return 0.0u"eV"
    end
    
    # 1. Generate all unique (i, j) pairs where i < j
    pairs = [(i, j) for i in 1:length(at) for j in (i+1):length(at)]
    
    # 2. Extract atomic positions as SVectors
    pos_cpu = [position(at, i) for i in 1:length(at)]
    
    # 3. Extract GPU-safe potential parameters
    pot_params = extract_gpu_params(pot)
    
    # 4. Unit cell matrix
    cell = SMatrix{3,3,Float64}(ustrip.(unitcell(at)))
    
    # 5. Transfer to GPU
    pairs_gpu = CuArray(pairs)
    pos_gpu = CuArray(pos_cpu)
    energy_gpu = CUDA.zeros(Float64, length(pairs))
    
    # 6. Launch kernel
    threads = 256
    blocks = cld(length(pairs), threads)
    CUDA.@sync @cuda threads=threads blocks=blocks compute_intra_pair_energies_kernel!(
        energy_gpu, pairs_gpu, pos_gpu, pot_params, cell
    )
    
    # 7. Return total energy with proper units
    return sum(Array(energy_gpu)) * u"eV"
end

# === GPU-Safe Helper Functions ===
@inline function pbc_dist_gpu(pos1::SVector{3, Float64}, pos2::SVector{3, Float64}, cell::SMatrix{3,3,Float64})
    dr = pos1 - pos2
    frac = inv(cell) * dr
    frac -= round.(frac)   # Wrap into unit cell
    return cell * frac     # Convert back to Cartesian
end

@inline function pair_energy_gpu(r::Float64, params::LJParams)
    # Safety check for very small distances
    if r < 1e-10
        return 0.0
    end
    σ = params.sigma
    ε = params.epsilon
    inv_r = σ / r
    inv_r6 = inv_r^6
    return 4.0 * ε * (inv_r6^2 - inv_r6)
end

# === Extract potential parameters in GPU format ===
function extract_gpu_params(pot::AbstractPotential)
    if pot isa LJParameters
        return LJParams(pot.epsilon.val, pot.sigma.val)
    else
        error("Unsupported potential type $(typeof(pot)). Please define `extract_gpu_params` for it.")
    end
end



"""
    frozen_energy(at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy of the frozen particles in the system using a composite Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `ljs::CompositeLJParameters{C}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, 
                       ljs::CompositeLJParameters{C}, 
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool}
                       ) where {C}
    check_num_components(C, list_num_par, frozen)
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(frozen) # find frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], ljs.lj_param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:C
        for j in 1:C
            if frozen[i] && frozen[j] && i != j # both frozen and different
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], ljs.lj_param_sets[i,j])
            end
        end
    end
    return energy
end

"""
    frozen_energy(at::AbstractSystem, lj::LJParameters, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy of the frozen particles in the system using a single Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, 
                       lj::LJParameters,
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool}
                       )
    if length(list_num_par) != length(frozen)
        throw(ArgumentError("The number of frozen and free parts does not match the length of the number of components."))
    end
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(frozen) # find frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], lj)
        end
    end
    # inter-component interactions
    for i in 1:length(list_num_par)
        for j in (i+1):length(list_num_par)
            if frozen[i] && frozen[j] # both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], lj)
            end
        end
    end
    return energy
end


"""
    interacting_energy(at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy from interactions between free-free and free-frozen particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `ljs::CompositeLJParameters{C}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy from interactions between particles.

"""
function interacting_energy(at::AbstractSystem, 
                            ljs::CompositeLJParameters{C}, 
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool}
                            ) where {C}
    check_num_components(C, list_num_par, frozen)
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], ljs.lj_param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:C
        for j in (i+1):C
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], ljs.lj_param_sets[i,j])
            end
        end
    end
    return energy
end

function interacting_energy(at::AbstractSystem, 
                            ljs::CompositeLJParameters{C}, 
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool},
                            surface::AbstractSystem
                            ) where {C}
    energy = 0.0u"eV"
    components_at = split_components(at, list_num_par)
    components = [components_at..., surface] # combine components from at and surface
    frozen = [frozen..., true] # add frozen state for surface
    list_num_par = [list_num_par..., length(surface)] # combine list_num_par vectors
    num_components = length(components)
    check_num_components(num_components, list_num_par, frozen)
    # @info "num_components: $num_components, list_num_par: $list_num_par, frozen: $frozen"
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], ljs.lj_param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:num_components
        for j in (i+1):num_components
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], ljs.lj_param_sets[i,j])
            end
        end
    end
    return energy
end

"""
    interacting_energy(at::AbstractSystem, lj::LJParameters, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy from interactions between free-free and free-frozen particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy from interactions between particles.

"""
function interacting_energy(at::AbstractSystem, 
                            lj::LJParameters,
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool}
                            )
    if length(list_num_par) != length(frozen)
        throw(ArgumentError("The number of frozen and free parts does not match the length of the number of components."))
    end
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], lj)
        end
    end
    # inter-component interactions
    for i in 1:length(list_num_par)
        for j in (i+1):length(list_num_par)
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], lj)
            end
        end
    end
    return energy
end

"""
    interacting_energy(at::AbstractSystem, lj::LJParameters)

Calculate the energy from interactions between particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.

# Returns
- `energy`: The energy from interactions between particles.

"""
interacting_energy(at::AbstractSystem, lj::LJParameters) = intra_component_energy(at, lj)

"""
    single_site_energy(index::Int, at::AbstractSystem, lj::LJParameters)
    single_site_energy(index::Int, at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int})
    single_site_energy(index::Int, at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int}, surface::AbstractSystem)

Calculate the energy of a single site in the system using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the site and all other sites in the system.

# Arguments
- `index::Int`: The index of the site for which the energy is calculated.
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `ljs::CompositeLJParameters{C}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `surface::AbstractSystem`: An optional surface system to consider in the energy calculation. See `LJSurfaceWalkers`.

# Returns
- `energy`: The energy of the site.

"""
function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            lj::LJParameters,
                            list_num_par::Vector{Int}
                            )
    
    all_index = collect(1:length(at))
    popat!(all_index, index)
    energies = Array{typeof(0.0u"eV"), 1}(undef, length(all_index))
    Threads.@threads for i in eachindex(all_index)
        r = pbc_dist(position(at, index), position(at, all_index[i]), at)
        energies[i] = lj_energy(r,lj)
    end
    return sum(energies)
end

function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            ljs::CompositeLJParameters{C},
                            list_num_par::Vector{Int},
                            ) where {C}
    comp_cut = vcat([0],cumsum(list_num_par))
    # @info "comp_cut: $comp_cut"
    comp_split = [comp_cut[i]+1:comp_cut[i+1] for i in 1:length(list_num_par)]
    # @info "comp_split: $comp_split"
    from_comp = findfirst(x->index in x, comp_split)
    # @info "from_comp: $from_comp"
    all_index = collect(1:length(at))
    popat!(all_index, index)
    energies = Array{typeof(0.0u"eV"), 1}(undef, length(all_index))
    Threads.@threads for i in eachindex(all_index)
        r = pbc_dist(position(at, index), position(at, all_index[i]), at)
        if all_index[i] in comp_split[from_comp]
            energy = lj_energy(r,ljs.lj_param_sets[from_comp,from_comp])
            energies[i] = energy
        else
            to_comp = findfirst(x->all_index[i] in x, comp_split)
            energy = lj_energy(r,ljs.lj_param_sets[from_comp,to_comp])
            energies[i] = energy
        end
    end
    return sum(energies)
end

function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            ljs::CompositeLJParameters{C},
                            list_num_par::Vector{Int},
                            surface::AbstractSystem
                            ) where {C}
    comp_cut = vcat([0],cumsum(list_num_par))
    # @info "comp_cut: $comp_cut"
    comp_split = [comp_cut[i]+1:comp_cut[i+1] for i in 1:length(list_num_par)]
    # @info "comp_split: $comp_split"
    from_comp = findfirst(x->index in x, comp_split)
    # @info "from_comp: $from_comp"
    all_index = collect(1:length(at))
    popat!(all_index, index)
    energies = Array{typeof(0.0u"eV"), 1}(undef, length(all_index))
    Threads.@threads for i in eachindex(all_index)
        r = pbc_dist(position(at, index), position(at, all_index[i]), at)
        if all_index[i] in comp_split[from_comp]
            energy = lj_energy(r,ljs.lj_param_sets[from_comp,from_comp])
            energies[i] = energy
        else
            to_comp = findfirst(x->all_index[i] in x, comp_split)
            energy = lj_energy(r,ljs.lj_param_sets[from_comp,to_comp])
            energies[i] = energy
        end
    end
    internal_e = sum(energies)
    energies_surface = Array{typeof(0.0u"eV"), 1}(undef, length(surface))
    Threads.@threads for i in eachindex(surface.position)
        r = pbc_dist(position(at, index), position(surface, i), at)
        to_comp = C # surface is the last component
        energies_surface[i] = lj_energy(r,ljs.lj_param_sets[from_comp,to_comp])
    end
    external_e = sum(energies_surface)
    return internal_e + external_e
end
