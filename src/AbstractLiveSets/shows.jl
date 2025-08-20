function Base.show(io::IO, walkers::LJAtomWalkers)
    println(io, "LJAtomWalkers($(eltype(walkers.walkers)), $(typeof(walkers.lj_potential))):")
    if length(walkers.walkers) > 10
        for i in 1:5
            println(io, "[$i] ", walkers.walkers[i])
        end
        println(io, "⋮\nOmitted ", length(walkers.walkers)-10, " walkers\n⋮\n")
        for i in length(walkers.walkers)-4:length(walkers.walkers)
            println(io, "[$i] ", walkers.walkers[i])
        end
    else
        for (ind, w) in enumerate(walkers.walkers)
            println(io, "[$ind] ", w)
        end
    end
    println(io, walkers.lj_potential)
end

function Base.show(io::IO, walkers::LJSurfaceWalkers)
    println(io, "LJSurfaceWalkers($(eltype(walkers.walkers)), $(typeof(walkers.potential))):")
    if length(walkers.walkers) > 10
        for i in 1:5
            println(io, "[$i] ", walkers.walkers[i])
        end
        println(io, "⋮\nOmitted ", length(walkers.walkers)-10, " walkers\n⋮\n")
        for i in length(walkers.walkers)-4:length(walkers.walkers)
            println(io, "[$i] ", walkers.walkers[i])
        end
    else
        for (ind, w) in enumerate(walkers.walkers)
            println(io, "[$ind] ", w)
        end
    end
    println(io, walkers.potential)
    println(io, "Surface: ", walkers.surface)
end