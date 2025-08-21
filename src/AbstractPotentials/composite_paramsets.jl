"""
    struct CompositeParameterSets{C,P} <: MultiComponentParameterSets

CompositeParameterSets is a struct that represents a set of composite parameter sets, typically used for multi-component systems.

# Fields
- `param_sets::Matrix{P}`: A matrix of parameter sets representing the composite parameter sets.

# Type Parameters
- `C::Int`: The number of composite parameter sets.

"""
struct CompositeParameterSets{C,P} <: MultiComponentParameterSets where {C,P}
    param_sets::Matrix{P}
    function CompositeParameterSets{C}(param_sets::Matrix{P}) where {C,P}
        if size(param_sets) != (C, C)
            throw(ArgumentError("The size of the matrix is not compatible with the number of components."))
        end
        new{C,P}(param_sets)
    end
end

function Base.show(io::IO, ps::CompositeParameterSets{C,P}) where {C,P}
    println(io, "CompositeParameterSets{$C}(param_sets::$(C)x$C Matrix{$P}):")
    for i in 1:C
        for j in 1:C
            println(io, "    param_sets[$i, $j] : ", ps.param_sets[i, j])
        end
    end
end

function Base.:(==)(cps1::CompositeParameterSets{C,P}, cps2::CompositeParameterSets{C,P}) where {C,P}
    return cps1.param_sets == cps2.param_sets
end

"""
    CompositeParameterSets(c::Int, ps::Vector{P}) where {P <: SingleComponentParameterSets}

Construct a `CompositeParameterSets` object from a vector of `SingleComponentParameterSets`.

# Arguments
- `c::Int`: The number of components.
- `ps::Vector{LJParameters}`: A vector of LJParameters. 
The number of elements in the vector must be equal to `c^2` or `c*(c+1)/2`. 
The former case is for a full flattened matrix of LJParameters, useful when 
the interactions are asymmetric, i.e., `epsilon_ij != epsilon_ji`. The latter
case is for symmetric interactions, i.e., `epsilon_ij = epsilon_ji`, hence only
the upper triangular part of the matrix is needed.

# Returns
- A `CompositeParameterSets` object.

# Example
```jldoctest
julia> ps = [LJParameters(epsilon=e) for e in [11, 12, 13, 22, 23, 33]]
6-element Vector{LJParameters}:
 LJParameters(11.0 eV, 1.0 Å, Inf, 0.0 eV)
 LJParameters(12.0 eV, 1.0 Å, Inf, 0.0 eV)
 LJParameters(13.0 eV, 1.0 Å, Inf, 0.0 eV)
 LJParameters(22.0 eV, 1.0 Å, Inf, 0.0 eV)
 LJParameters(23.0 eV, 1.0 Å, Inf, 0.0 eV)
 LJParameters(33.0 eV, 1.0 Å, Inf, 0.0 eV)

julia> ljp = CompositeParameterSets(3,ps)
CompositeParameterSets{3}(param_sets::3x3 Matrix{LJParameters}):
    param_sets[1, 1] : LJParameters(11.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[1, 2] : LJParameters(12.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[1, 3] : LJParameters(13.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[2, 1] : LJParameters(12.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[2, 2] : LJParameters(22.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[2, 3] : LJParameters(23.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[3, 1] : LJParameters(13.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[3, 2] : LJParameters(23.0 eV, 1.0 Å, Inf, 0.0 eV)
    param_sets[3, 3] : LJParameters(33.0 eV, 1.0 Å, Inf, 0.0 eV)
```
"""
function CompositeParameterSets(c::Int, ps::Vector{P}) where {P <: SingleComponentParameterSets}
    if length(ps) == c^2
        # If the number of LJParameters sets is equal to the number 
        # of elements in the matrix, then assume that Vector{LJParameters} 
        # is a flattened matrix and reshape the vector into a matrix.
        @info "Creating CompositeParameterSets from a flattened matrix. By specifying
        $length(ps) sets of parameters, a $c x $c matrix is constructed. If this
        was not your intention, please check the documentation or raise an issue."
        return CompositeParameterSets{c}(reshape(ps, c, c))
    elseif length(ps) == c*(c+1)/2
        # If the number of LJParameters sets is equal to the number
        # of elements in the upper triangular part of the matrix, then
        # construct the full matrix from the upper triangular part.
        @info "Creating CompositeParameterSets from the upper triangular part of the matrix.
        By specifying $length(ps) sets of parameters, a $c x $c matrix is constructed.
        If this was not your intention, please check the documentation or raise an issue."
        ljmatrix = Matrix{P}(undef, c, c)
        k = 1
        for i in 1:c
            for j in i:c
            ljmatrix[i, j] = ps[k]
            if i != j
                ljmatrix[j, i] = ps[k]
            end
            k += 1
            end
        end
        return CompositeParameterSets{c}(ljmatrix)
    else
        throw(ArgumentError("the number of parameter sets is not compatible with the number of components."))
    end
end