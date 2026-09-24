# The random stream of the grand-canonical nested-sampling kernels is independent of the
# energies: the number of random draws a kernel step consumes never depends on the ceiling
# outcome and never on the acceptance ratio. Every insertion and deletion proposal draws its
# Metropolis uniform immediately after its proposal draws, whatever the ceiling decides and
# whatever the ratio, so a walk that draws its trial proposals before evaluating any of
# their energies consumes exactly the stream of the serial kernel. The tests below take
# one step from a seeded generator under a ceiling that accepts everything and under one
# that rejects everything, at an acceptance ratio far above one and far below one, and
# read the next number of the stream against a probe of the same seed: the index of the
# probe entry it equals is the number of draws the step consumed. The particle count after
# the step is asserted beside the draw count, so that every arm is shown to reach the branch
# it names: an accepted insertion at the high activity, a rejected one at the low activity,
# and the reverse for a deletion.
#
# Channel selection reads one uniform r: r < p_move is a displacement (a swap on the
# lattice), p_move <= r < p_move + p_insert an insertion, the rest a deletion. The mixes
# (p_move, p_insert) = (0, 0.75) and (0, 0.25) keep the insertion and the deletion channel
# both open, so the acceptance ratios z0V (p_delete / p_insert) / (n + 1) and
# (p_insert / p_delete) n / z0V (their lattice forms carry the same factors) are nonzero and
# the activities 1e6 and 1e-6 put them far above and far below one; a seed whose first draw
# lies in [0.25, 0.75) selects the insertion under the first mix and the deletion under the
# second, which is asserted as a precondition. A closed reverse channel (p_delete = 0 or
# p_insert = 0) makes every ratio exactly zero: it forces a rejection but never reaches the
# ratio-at-or-above-one branch, where the conditional draw used to draw nothing.
#
# Zero-interaction fixtures (epsilon = 0): every energy is exactly 0.0 eV, so the ceilings
# +1 eV and -1 eV decide every proposal the same way on every Julia version and
# architecture; the counts are exact integers.
@testset "grand-canonical kernels: the draw count of a step never depends on an energy" begin
    using Random

    # the probe: the first entries of the stream at this seed; a step that consumed k draws
    # leaves the (k + 1)-th entry as the next number
    probe(seed, n=16) = (Random.seed!(seed); [rand() for _ in 1:n])
    function consumed(seed, f, n=16)
        p = probe(seed, n)
        Random.seed!(seed)
        f()
        x = rand()
        k = findfirst(==(x), p)
        return k === nothing ? -1 : k - 1
    end
    ceilings = (1.0u"eV", -1.0u"eV")            # accept everything, reject everything
    activities = (1.0e6, 1.0e-6)                # ratio far above one, far below one
    seeds = (778, 4243)                         # first draws 0.3139 and 0.5068
    for seed in seeds
        @test 0.25 <= probe(seed, 1)[1] < 0.75  # branch precondition: insertion under (0, 0.75), deletion under (0, 0.25)
    end
    # the particle count after one step from n = 2: an insertion is accepted only at the high
    # activity, a deletion only at the low one, and nothing changes under the rejecting ceiling
    function n_after(label, accept::Bool, z0V)
        accept || return 2
        label == "insertion" && return z0V > 1.0 ? 3 : 2
        label == "deletion" && return z0V > 1.0 ? 2 : 1
        return 2
    end

    @testset "atomistic kernel, pairwise potential" begin
        box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, (true, true, true)))
        mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                               empty(position(seed_at, :)), empty(species(seed_at, :)),
                               empty(mass(seed_at, :)))
        lj0 = LJParameters(epsilon=0.0, sigma=2.5, cutoff=2.5)
        function mkwalker(n)
            w = AtomWalker{1}(mkempty())
            for i in 1:n
                insert_particle!(w, SVector(2.0 * i, 5.0, 5.0)u"Å", :Ar)
            end
            return w
        end
        # (p_move, p_insert) = (0, 0.75) selects an insertion and (0, 0.25) a deletion at these
        # seeds; (1, 0) forces a displacement
        for (label, mix, expected) in (("insertion", (0.0, 0.75), 5),     # channel, x, y, z, uniform
                                       ("deletion", (0.0, 0.25), 3),      # channel, index, uniform
                                       ("displacement", (1.0, 0.0), 5))   # channel, index, dx, dy, dz
            for emax in ceilings, z0V in activities, seed in seeds
                w = mkwalker(2)
                k = consumed(seed, () -> MC_grand_canonical_walk!(1, w, lj0, emax; z0V=z0V, species=:Ar,
                                                                  p_move=mix[1], p_insert=mix[2]))
                @test k == expected
                @test w.list_num_par == [n_after(label, emax > 0.0u"eV", z0V)]
            end
        end
        # guard skips consume the channel draw only, at both ceilings
        for emax in ceilings
            w = mkwalker(0)
            @test consumed(777, () -> MC_grand_canonical_walk!(1, w, lj0, emax; z0V=8.0, species=:Ar,
                                                              p_move=0.0, p_insert=0.0)) == 1
            w = mkwalker(0)
            @test consumed(777, () -> MC_grand_canonical_walk!(1, w, lj0, emax; z0V=8.0, species=:Ar,
                                                              p_move=1.0, p_insert=0.0)) == 1
            w = mkwalker(1)
            @test consumed(777, () -> MC_grand_canonical_walk!(1, w, lj0, emax; z0V=8.0, species=:Ar,
                                                              p_move=0.0, p_insert=1.0, n_max=1)) == 1
        end
        # the cavity-biased channel: the sub-channel draw, then the cell draw and three
        # jitters (biased) or the three uniforms, then the Metropolis uniform, at both
        # ceilings and both activities
        for emax in ceilings, z0V in activities, seed in seeds
            w = mkwalker(2)
            kb = consumed(seed, () -> MC_grand_canonical_walk!(1, w, lj0, emax; z0V=z0V, species=:Ar,
                                                               p_move=0.0, p_insert=0.75, p_bias=1.0,
                                                               bias_radius=1.0, bias_grid=4))
            @test kb == 7                          # channel, sub-channel, cell, three jitters, uniform
            @test w.list_num_par == [n_after("insertion", emax > 0.0u"eV", z0V)]
            w = mkwalker(2)
            ku = consumed(seed, () -> MC_grand_canonical_walk!(1, w, lj0, emax; z0V=z0V, species=:Ar,
                                                               p_move=0.0, p_insert=0.75, p_bias=1.0e-12,
                                                               bias_radius=1.0, bias_grid=4))
            @test ku == 6                          # channel, sub-channel, x, y, z, uniform
            @test w.list_num_par == [n_after("insertion", emax > 0.0u"eV", z0V)]
        end
        # a walk of many steps: the stream position after the walk is the same under both
        # ceilings whenever the channel sequence is (deletions at N = 0 excluded by starting
        # with the particles a rejected walk keeps): fixed-N mix of displacements only
        for seed in (4241, 4242)
            positions = Int[]
            for emax in ceilings
                w = mkwalker(3)
                push!(positions, consumed(seed, () -> MC_grand_canonical_walk!(11, w, lj0, emax; z0V=3.0,
                                                                                species=:Ar, p_move=1.0, p_insert=0.0), 64))
            end
            @test positions[1] == positions[2] == 11 * 5
        end
    end

    @testset "atomistic kernel, surface-aware method" begin
        sbox = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
        surf_sys = FastSystem(atomic_system([:H => [2.5, 2.5, 2.0]u"Å", :H => [7.5, 2.5, 2.0]u"Å",
                                             :H => [2.5, 7.5, 2.0]u"Å", :H => [7.5, 7.5, 2.0]u"Å"],
                                            sbox, (true, true, false)))
        surf = AtomWalker(deepcopy(surf_sys); freeze_species=[:H])
        smkempty() = FastSystem(cell_vectors(surf_sys), periodicity(surf_sys),
                                empty(position(surf_sys, :)), empty(species(surf_sys, :)),
                                empty(mass(surf_sys, :)))
        lj0 = LJParameters(epsilon=0.0, sigma=2.5, cutoff=1.8, shift=true)
        cps0 = CompositeParameterSets(2, [lj0, lj0, lj0])
        function smkwalker(n)
            w = AtomWalker{1}(smkempty())
            for i in 1:n
                insert_particle!(w, SVector(2.0 * i, 5.0, 8.0)u"Å", :Ar)
            end
            return w
        end
        for (label, mix, expected) in (("insertion", (0.0, 0.75), 5), ("deletion", (0.0, 0.25), 3),
                                       ("displacement", (1.0, 0.0), 5))
            for emax in ceilings, z0V in activities, seed in seeds
                w = smkwalker(2)
                k = consumed(seed, () -> MC_grand_canonical_walk!(1, w, cps0, emax, surf; z0V=z0V,
                                                                  species=:Ar, p_move=mix[1], p_insert=mix[2]))
                @test k == expected
                @test w.list_num_par == [n_after(label, emax > 0.0u"eV", z0V)]
            end
        end
        for emax in ceilings
            w = smkwalker(0)
            @test consumed(777, () -> MC_grand_canonical_walk!(1, w, cps0, emax, surf; z0V=8.0, species=:Ar,
                                                              p_move=0.0, p_insert=0.0)) == 1
        end
    end

    @testset "lattice kernel" begin
        # a 4 x 4 square lattice with two occupied sites; the zero Hamiltonian makes every
        # energy exactly 0.0 eV, the perturbation is left at its default of zero
        occ = falses(16); occ[1] = true; occ[6] = true
        lat = SLattice{SquareLattice}(supercell_dimensions=(4, 4, 1), components=[collect(occ)])
        h0 = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        omegas = (1.0, -1.0)
        occupied(lw) = sum(lw.configuration.components[1])
        for (label, mix, expected) in (("insertion", (0.0, 0.75), 4),   # channel, site, uniform, perturbation
                                       ("deletion", (0.0, 0.25), 4),    # channel, site, uniform, perturbation
                                       ("swap", (1.0, 0.0), 4))         # channel, two sites, perturbation
            for om in omegas, z0 in activities, seed in seeds
                lw = LatticeWalker(deepcopy(lat))
                k = consumed(seed, () -> MC_grand_canonical_walk!(1, lw, h0, om, 0.0; p_move=mix[1],
                                                                  p_insert=mix[2], z0=z0))
                @test k == expected
                @test occupied(lw) == n_after(label, om > 0.0, z0)
            end
        end
        # the biased insertion channel: channel, sub-channel, site, uniform, perturbation
        for om in omegas, z0 in activities, seed in seeds
            lw = LatticeWalker(deepcopy(lat))
            k = consumed(seed, () -> MC_grand_canonical_walk!(1, lw, h0, om, 0.0; p_move=0.0, p_insert=0.75,
                                                              z0=z0, p_bias=0.5, bias_predicate=:contact))
            @test k == 5
            @test occupied(lw) == n_after("insertion", om > 0.0, z0)
        end
        # guard skips: an insertion at the cap and a deletion on the empty lattice consume
        # the channel draw only
        for om in omegas
            lw = LatticeWalker(deepcopy(lat))
            @test consumed(4243, () -> MC_grand_canonical_walk!(1, lw, h0, om, 0.0; p_move=0.0, p_insert=1.0,
                                                               n_max=2)) == 1
            empty_lat = SLattice{SquareLattice}(supercell_dimensions=(4, 4, 1), components=[collect(falses(16))])
            lw = LatticeWalker(empty_lat)
            @test consumed(4242, () -> MC_grand_canonical_walk!(1, lw, h0, om, 0.0; p_move=0.0, p_insert=0.0)) == 1
        end
    end
end
