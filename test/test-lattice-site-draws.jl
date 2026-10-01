# The site draws of the grand-canonical lattice moves (the empty site of an insertion, the
# occupied site of a deletion, the two sites of an occupied-to-empty hop, and the inlined
# uniform draws of the biased insertion channel) go through `_rand_site(occupancy, occupied,
# count)`: one `rand(1:count)` draw picks a rank and a scan returns the site of that rank in
# index order, where the moves used to draw `rand(findall(...))`. `rand(v::Vector)` indexes
# `v` with the same `rand(1:length(v))`, so the draw, the site and the rest of the random
# stream are unchanged. The seeded grand-canonical tests elsewhere in the suite hold the
# number of draws of every step; which site each call draws is checked here: draw by draw
# for the helper against the expression it replaces, on occupancy vectors of 1, 2, 5, 64 and
# 256 sites at every filling class (empty, one site, half, all but one, full), as
# `Vector{Bool}` (a lattice's occupancy) and as `BitVector` (the helper takes any
# `AbstractVector{Bool}`); step by step for the exported insertion and deletion; and move by
# move for the occupied-to-empty hop of `MC_grand_canonical_walk!`, against a replay of the
# walk with the replaced expressions. The biased channel's inlined draws call the same
# helper with the same counts as the exported functions. Then the allocation: the bytes per
# move of a walk with occupied-to-empty hops, insertions, deletions and site deltas on a
# 32 x 32 lattice are at most 512 above those on an 8 x 8 lattice (they grew from about 0.7
# to 9.3 kB per move, one index vector of up to M entries per draw), and a single draw
# allocates nothing. The biased channel still allocates in proportion to the lattice through
# `lattice_biased_sites`, so its walk is not part of the allocation check.
@testset "lattice site draws: the draws of rand(findall(...)), without allocating" begin
    using Random

    @testset "one draw against rand(findall(...)), stream position included" begin
        for M in (1, 2, 5, 64, 256), container in (Vector{Bool}, BitVector)
            for n in sort(unique([0, 1, M ÷ 2, M - 1, M]))
                Random.seed!(1000 * M + n)
                occupancy = container(fill(false, M))
                occupancy[randperm(M)[1:n]] .= true
                for occupied in (true, false)
                    count = occupied ? n : M - n
                    count == 0 && continue
                    for seed in 1:8
                        Random.seed!(seed)
                        expected = rand(findall(occupancy .== occupied))
                        expected_next = rand()
                        Random.seed!(seed)
                        @test FreeBird.MonteCarloMoves._rand_site(occupancy, occupied, count) == expected
                        @test rand() == expected_next
                    end
                end
            end
        end
    end

    @testset "lattice_insert_particle! and lattice_delete_particle! step by step" begin
        # alternate insertions and deletions from an empty, a half-filled and a full lattice;
        # the reference replays the same steps with the expressions the functions used to draw
        for L in (4, 8), filling in (0, 1, 2), seed in 1:4
            Random.seed!(L + 10 * filling)
            occupancy = fill(false, L^2)
            occupancy[randperm(L^2)[1:(filling * L^2 ÷ 2)]] .= true
            lattice = SLattice{SquareLattice}(supercell_dimensions=(L, L, 1), components=[copy(occupancy)])
            Random.seed!(seed)
            sites = [isodd(call) ? lattice_insert_particle!(lattice)[3] : lattice_delete_particle!(lattice)[3]
                     for call in 1:40]
            next = rand()
            Random.seed!(seed)
            expected = Int[]
            for call in 1:40
                if isodd(call)
                    site = all(occupancy) ? 0 : rand(findall(.!occupancy))
                    site > 0 && (occupancy[site] = true)
                else
                    site = any(occupancy) ? rand(findall(occupancy)) : 0
                    site > 0 && (occupancy[site] = false)
                end
                push!(expected, site)
            end
            @test sites == expected
            @test lattice.components[1] == occupancy
            @test rand() == next
        end
    end

    @testset "the occupied-to-empty hop of MC_grand_canonical_walk! move by move" begin
        # hops only (p_move = 1), an open ceiling: every step draws the channel uniform, then
        # (unless the lattice is empty or full, a guard skip) one occupied and one empty site,
        # exchanges them, and draws the tie-breaking uniform; the replay draws the same numbers
        # with the replaced expressions
        h = GenericLatticeHamiltonian(0.0, [-0.05, 0.0], u"eV")
        function replay!(occupancy, steps)
            attempted = 0
            for _ in 1:steps
                rand()                                   # the channel draw
                n = count(occupancy)
                (n == 0 || n == length(occupancy)) && continue
                hop_from = rand(findall(occupancy))
                hop_to = rand(findall(.!occupancy))
                occupancy[hop_from], occupancy[hop_to] = false, true
                attempted += 1
                rand()                                   # the tie-breaking draw
            end
            return attempted
        end
        for L in (4, 8, 16), filling in (0, 1, 2, 3, 4), seed in 1:3, incremental in (false, true)
            M = L^2
            n = (0, 1, M ÷ 2, M - 1, M)[filling+1]
            Random.seed!(100 * L + seed)
            start = fill(false, M)
            start[randperm(M)[1:n]] .= true
            walker = LatticeWalker(SLattice{SquareLattice}(supercell_dimensions=(L, L, 1), components=[copy(start)]))
            assign_energy!(walker, h)
            Random.seed!(seed)
            counters = MC_grand_canonical_walk!(300, walker, h, Inf, 0.0; p_move=1.0, p_insert=0.0,
                                                swap_mode=:occupied_empty, incremental=incremental)[end]
            next = rand()
            occupancy = copy(start)
            Random.seed!(seed)
            attempted = replay!(occupancy, 300)
            @test walker.configuration.components[1] == occupancy
            @test rand() == next
            @test counters.swap_attempted == counters.swap_accepted == attempted
        end
    end

    @testset "allocation that does not grow with the lattice" begin
        h = GenericLatticeHamiltonian(0.0, [-0.05, 0.0], u"eV")
        bytes_per_move = Dict{Int,Float64}()
        for L in (8, 32)
            Random.seed!(L)
            occupancy = fill(false, L^2)
            occupancy[randperm(L^2)[1:(L^2 ÷ 2)]] .= true
            walker = LatticeWalker(SLattice{SquareLattice}(supercell_dimensions=(L, L, 1), components=[occupancy]))
            assign_energy!(walker, h)
            walk!(n) = MC_grand_canonical_walk!(n, walker, h, Inf, 0.0; p_move=0.4, p_insert=0.3,
                                                swap_mode=:occupied_empty, incremental=true)
            # a long first walk reaches every channel and many rejections, so every method the
            # walk dispatches at run time is compiled before the measured walks
            walk!(4000)
            # the slope between walks of 2,000 and 4,000 moves: the fixed cost of a call cancels
            short = @allocated walk!(2000)
            long = @allocated walk!(4000)
            bytes_per_move[L] = (long - short) / 2000
        end
        @test bytes_per_move[32] <= bytes_per_move[8] + 512

        occupancy = fill(false, 256)
        occupancy[1:2:end] .= true
        draw(occupancy) = FreeBird.MonteCarloMoves._rand_site(occupancy, false, 128)
        draw(occupancy)
        @test (@allocated draw(occupancy)) == 0
    end

    @testset "counts the occupancy cannot meet" begin
        # a count of zero is an empty range, as findall(...) of no site was an empty vector
        @test_throws ArgumentError FreeBird.MonteCarloMoves._rand_site([false, true], true, 0)
        # with no matching site, any count runs off the end of the scan
        @test_throws ArgumentError FreeBird.MonteCarloMoves._rand_site([false, false], true, 1)
    end
end
