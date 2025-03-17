"""
    periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T

Wrap the position vector `pos` according to the periodic boundary conditions of the `system`. If the
boundary condition is `Periodic()`, the position is wrapped using the modulo operator. If the boundary
condition is `DirichletZero()`, the position is wrapped by reflecting the position vector across the boundary.

# Arguments
- `pos::SVector{3,T}`: The position vector to be wrapped.
- `system::AbstractSystem`: The system containing the periodic boundary conditions.

# Returns
The wrapped position vector.

"""
function periodic_boundary_wrap!(pos::SVector{3,T}, system::AbstractSystem) where T
    pbc = periodicity(system)
    box = cell_vectors(system)
    new_pos = Vector{typeof(0.0u"Å")}(undef, 3)
    for i in eachindex(pbc)
        if pbc[i] == true # wrap the position
            new_pos[i] = mod(pos[i], box[i][i])
        elseif pbc[i] == false # reflect the position
            if pos[i] > box[i][i]
                new_pos[i] = box[i][i]*2 - pos[i]
            elseif pos[i] < 0.0u"Å"
                new_pos[i] = -pos[i]
            else
                new_pos[i] = pos[i]
            end
        else
            error("Unsupported boundary condition: $(pbc[i]).")
        end
    end
    pos = SVector{3,T}(new_pos)
    return pos
end

### EXPERIMENTAL FUNCTION! ###
"""
    mean_sq_displacement(at::AtomWalker, at_orig::AtomWalker)

Calculate the mean squared displacement before and after random walk(s). Note that due to the current implementation of the periodic boundary wrap, this function is not appropriate to use for calculating mean displacements in a propagation.

# Arguments
- `at::AtomWalker{C}`: The current `AtomWalker` after the random walk.
- `at_orig::AtomWalker{C}`: The original `AtomWalker` before the random walk.

# Returns
- `distsq::typeof(0.0u"Å"^2)`: The mean squared displacement of all free particles.

"""
function mean_sq_displacement(at::AtomWalker{C}, at_orig::AtomWalker{C}) where C
    distsq::typeof(0.0u"Å"^2) = 0.0u"Å"^2
    free_index = free_par_index(at)
    for i in free_index
        dist::typeof(0.0u"Å") = pbc_dist(at.configuration.position[i], at_orig.configuration.position[i], at.configuration)
        distsq += dist^2
    end
    return distsq/length(free_index)
end

"""
    free_par_index(at::AtomWalker{C}) where C

Get the indices of the free particles in the `AtomWalker`.
"""
function free_par_index(at::AtomWalker{C}) where C
    ind_free_par = Array{Int}(undef, 0)
    comp_cut = vcat([0],cumsum(at.list_num_par))
    # @show comp_cut
    for i in 1:C
        if !at.frozen[i]
            ind_free_par = vcat(ind_free_par,comp_cut[i]+1:comp_cut[i+1])
        end
    end
    return ind_free_par
end


"""
    free_component_index(at::AtomWalker{C}) where C

Get the indices of the free particles in each component of the `AtomWalker`.

# Returns
- `ind_free_parts::Array{Vector{Int}}`: An array of vectors containing the indices of the free particles in each component.
"""
function free_component_index(at::AtomWalker{C}) where C
    ind_free_parts = Array{Vector{Int}}(undef, 0)
    comp_cut = vcat([0],cumsum(at.list_num_par))
    # @show comp_cut
    for i in 1:C
        if !at.frozen[i]
            push!(ind_free_parts,collect(comp_cut[i]+1:comp_cut[i+1]))
        end
    end
    return ind_free_parts
end