@testset "Pairwise forces and displacement vectors" begin
    using Random

    box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    lj = LJParameters(epsilon=0.05, sigma=2.5, cutoff=2.0, shift=true)

    @testset "pair_force: analytic vs finite difference" begin
        # well, wall, and mid-range; the FD fallback and the analytic override
        # agree to the FD truncation class away from the cutoff
        for r in (2.2, 2.5, 2.806, 3.2, 4.0, 4.8)
            f_fd = ustrip(u"eV/Å",
                invoke(pair_force, Tuple{typeof(1.0u"Å"), SingleComponentPotential{Pairwise}},
                       (r)u"Å", lj))
            f_an = ustrip(u"eV/Å", pair_force((r)u"Å", lj))
            # rtol for the generic scale; atol covers the zero-crossing
            # neighborhood, where the FD absolute error class is ~1e-10
            @test isapprox(f_an, f_fd; rtol=1e-7, atol=1e-9)
        end
        # the truncation discontinuity is pinned from both sides
        rc = 2.0 * 2.5
        @test pair_force((rc + 1e-9)u"Å", lj) == 0.0u"eV/Å"
        @test abs(ustrip(u"eV/Å", pair_force((rc - 1e-9)u"Å", lj))) > 1e-4
        # zero crossing at the minimum r0 = 2^(1/6) sigma
        r0 = 2.5 * 2.0^(1 / 6)
        @test abs(ustrip(u"eV/Å", pair_force((r0)u"Å", lj))) < 1e-12
        @test ustrip(u"eV/Å", pair_force((r0 - 0.1)u"Å", lj)) > 0.0
        @test ustrip(u"eV/Å", pair_force((r0 + 0.1)u"Å", lj)) < 0.0
    end

    @testset "pbc_displacement: minimum image and norm agreement" begin
        pbc_sys = FastSystem(atomic_system([:Ar => [0.5, 5.0, 5.0]u"Å",
                                            :Ar => [9.5, 5.0, 5.0]u"Å"], box, (true, true, true)))
        d = pbc_displacement(position(pbc_sys, 1), position(pbc_sys, 2), pbc_sys)
        @test isapprox(ustrip(u"Å", d[1]), 1.0; atol=1e-12)   # wraps, not -9
        @test ustrip(u"Å", d[2]) == 0.0 && ustrip(u"Å", d[3]) == 0.0
        mixed = FastSystem(atomic_system([:Ar => [0.5, 5.0, 0.5]u"Å",
                                          :Ar => [9.5, 5.0, 9.5]u"Å"], box, (true, true, false)))
        dm = pbc_displacement(position(mixed, 1), position(mixed, 2), mixed)
        @test isapprox(ustrip(u"Å", dm[1]), 1.0; atol=1e-12)   # periodic axis wraps
        @test isapprox(ustrip(u"Å", dm[3]), -9.0; atol=1e-12)  # open axis does not
        Random.seed!(30301)
        for _ in 1:50
            sys = FastSystem(periodic_system([:Ar => rand(3), :Ar => rand(3)], box, fractional=true))
            dd = pbc_displacement(position(sys, 1), position(sys, 2), sys)
            nrm = sqrt(sum(x -> ustrip(u"Å", x)^2, dd))
            @test isapprox(nrm, ustrip(u"Å", pbc_dist(position(sys, 1), position(sys, 2), sys));
                           rtol=1e-12, atol=1e-12)
        end
    end

    @testset "interacting_gradient vs finite differences" begin
        Random.seed!(30302)
        coor = [:Ar => rand(3) for _ in 1:5]
        sys = FastSystem(periodic_system(coor, box, fractional=true))
        w = AtomWalker(sys)
        g = interacting_gradient(w.configuration, lj, w.list_num_par, w.frozen)
        h = 1e-5
        for i in 1:5, k in 1:3
            orig = position(w.configuration, i)
            shift = SVector(k == 1 ? h : 0.0, k == 2 ? h : 0.0, k == 3 ? h : 0.0)u"Å"
            w.configuration.position[i] = orig + shift
            Ep = ustrip(u"eV", interacting_energy(w.configuration, lj, w.list_num_par, w.frozen))
            w.configuration.position[i] = orig - shift
            Em = ustrip(u"eV", interacting_energy(w.configuration, lj, w.list_num_par, w.frozen))
            w.configuration.position[i] = orig
            @test isapprox(ustrip(u"eV/Å", g[i][k]), (Ep - Em) / (2h); rtol=1e-6, atol=1e-10)
        end
        # translation invariance: the gradient sums to zero over an all-free
        # periodic system
        tot = sum(g)
        @test all(abs(ustrip(u"eV/Å", tot[k])) < 1e-12 for k in 1:3)
        # symmetric fixture: the middle particle of an equidistant colinear
        # triple feels no net force
        trip = AtomWalker(FastSystem(atomic_system(
            [:Ar => [2.0, 5.0, 5.0]u"Å", :Ar => [4.8, 5.0, 5.0]u"Å",
             :Ar => [7.6, 5.0, 5.0]u"Å"], box, (true, true, true))))
        gt = interacting_gradient(trip.configuration, lj, trip.list_num_par, trip.frozen)
        @test all(abs(ustrip(u"eV/Å", gt[2][k])) < 1e-14 for k in 1:3)
    end

    @testset "frozen split" begin
        # one free Ar and two frozen Xe: the frozen entries are exactly zero,
        # the free entry matches the hand-built two-body sum, and the
        # frozen-frozen pair contributes nothing anywhere
        sys = FastSystem(atomic_system([:Ar => [3.0, 5.0, 5.0]u"Å",
                                        :Xe => [6.0, 5.0, 5.0]u"Å",
                                        :Xe => [6.0, 7.5, 5.0]u"Å"], box, (true, true, true)))
        w = AtomWalker(sys; freeze_species=[:Xe])
        @test w.frozen == [true, false] || w.frozen == [false, true]
        g = interacting_gradient(w.configuration, lj, w.list_num_par, w.frozen)
        # locate the free particle (components are sorted by atomic number)
        free_idx = findfirst(i -> species(w.configuration, i) == ChemicalSpecies(:Ar), 1:3)
        for i in 1:3
            if i == free_idx
                hand = SVector(0.0u"eV/Å", 0.0u"eV/Å", 0.0u"eV/Å")
                for j in 1:3
                    j == free_idx && continue
                    disp = pbc_displacement(position(w.configuration, free_idx),
                                            position(w.configuration, j), w.configuration)
                    r = sqrt(sum(x -> x^2, disp))
                    hand += (-pair_force(r, lj) / r) * disp
                end
                # atol per the house tie-charge convention; equally sharp for
                # a two-term sum, robust to cross-architecture reassociation
                @test all(isapprox(g[free_idx][k], hand[k]; atol=1e-14u"eV/Å") for k in 1:3)
            else
                @test g[i] == SVector(0.0u"eV/Å", 0.0u"eV/Å", 0.0u"eV/Å")
            end
        end
    end

    @testset "short-circuits" begin
        empty_sys = FastSystem(box, (true, true, true),
                               SVector{3, typeof(1.0u"Å")}[], ChemicalSpecies[],
                               typeof(1.0u"u")[])
        @test interacting_gradient(empty_sys, lj, [0], [false]) == SVector{3, typeof(0.0u"eV/Å")}[]
        single = FastSystem(atomic_system([:Ar => [5.0, 5.0, 5.0]u"Å"], box, (true, true, true)))
        g1 = interacting_gradient(single, lj, [1], [false])
        @test length(g1) == 1
        @test g1[1] == SVector(0.0u"eV/Å", 0.0u"eV/Å", 0.0u"eV/Å")
    end
end
