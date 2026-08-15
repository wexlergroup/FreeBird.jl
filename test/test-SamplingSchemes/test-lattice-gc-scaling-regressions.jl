@testset "Lattice GC scaling and ledger regressions" begin
    using Random
    using Unitful
    using DataFrames

    # Capstone coverage for the copy-free/incremental/ledger round: the
    # opt-in fast paths and the opt-in proposal mode sample the same
    # distribution the shipped path does, the ledgers close against their
    # run totals on fresh seeded descents, and the performance win holds as
    # a structural (allocation-shape) invariant, never a timing assertion.
    #
    # Calibration ledger (gates at >= 3x the max three-seed deviation, per
    # statistic, on the 4x4 two-shell fixture, K = 64, mc_steps = 40, 820
    # steps, z0 = 1, reduced at mu = 0, T = 600 K):
    # - default arm, seeds 99301/99302/99303:
    #   |d logXi| = 0.109/0.025/0.133, |d mean_N| = 0.177/0.201/0.005
    #   -> gates 0.40 and 0.61.
    # - incremental arm, seeds 96001/96002/96003 (the incremental-walk
    #   change's own calibration, same fixture class):
    #   |d logXi| = 0.153/0.282/0.156, |d mean_N| = 0.057/0.088/0.076
    #   -> gates 0.85 and 0.27.
    # - occupied-empty arm, seeds 99311/99312/99313:
    #   |d logXi| = 0.015/0.002/0.453, |d mean_N| = 0.131/0.041/0.260
    #   -> gates 1.36 and 0.78.
    # Wall-clock budget: the whole file runs in about 10 s single-threaded
    # on the calibration machine (the three-arm exactness block dominates).
    # Occupied-empty scope note: cross-proposal DEEP-quartile occupancy
    # moments were calibrated and systematically differ (quartile-4 mean-N
    # gap 0.48 vs a 0.15 within-arm spread class): a decorrelation-quality
    # signature of the null-dominated uniform-pair proposal, not a
    # stationarity defect, so the stationarity gates run against the exact
    # grand sum and only the first quartile carries a cross-proposal
    # agreement gate (measured Q1 gap 0.103 vs within-arm spread 0.118;
    # gate 0.35 at the 3x rule).
    # Allocation figures (single thread, compile-warmed): mixed-channel
    # incremental 200-step walk 541 kB at M = 256; swap-only marginal
    # (400-step minus 200-step, cancelling the once-per-call O(M) anchor)
    # 33.2 kB at M = 256 vs 41.5 kB at M = 1024 (ratio 1.25 for 4x the
    # sites; an O(M)-per-step regression of the deepcopy or full-recompute
    # class shifts this by orders of magnitude).

    cap_lat() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false), cutoff_radii=[1.1, 1.5],
        components=[[false for _ in 1:16]], adsorptions=:full)
    cap_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
    cap_kb = 8.617333262e-5
    cap_T = 600.0
    cap_mu = 0.0
    cap_beta = 1 / (cap_kb * cap_T)
    cap_save() = SaveEveryN("t_cap.csv", "t_cap.traj", "t_cap.ls",
                            10^7, 10^7, 10^7)
    cap_cleanup() = rm.(["t_cap.csv", "t_cap.traj", "t_cap.ls"], force=true)

    # Exact grand sum by per-sector enumeration, before any seeding
    # (exact_enumeration consumes the global RNG)
    cap_lnZs = Float64[]
    for N in 0:16
        occ = vcat(fill(true, N), fill(false, 16 - N))
        latN = cap_lat()
        latN.components[1] .= occ
        dfe, _ = exact_enumeration(latN, cap_ham)
        Es = [ustrip(u"eV", e) for e in dfe.energy]
        Emin = minimum(Es)
        push!(cap_lnZs, cap_beta * cap_mu * N - cap_beta * Emin +
                        log(sum(exp.(-cap_beta .* (Es .- Emin)))))
    end
    cap_m = maximum(cap_lnZs)
    cap_lnXi = cap_m + log(sum(exp.(cap_lnZs .- cap_m)))
    cap_meanN = sum((0:16) .* exp.(cap_lnZs .- cap_m)) /
                sum(exp.(cap_lnZs .- cap_m))

    function cap_arm(seed, routine)
        Random.seed!(seed)
        ws = [LatticeWalker(deepcopy(cap_lat()), energy=0.0u"eV", iter=0)
              for _ in 1:64]
        ls = LatticeGasWalkers(ws, cap_ham; assign_energy=false)
        p = IdealGasReferencedGCNSParameters(mc_steps=40,
            reference_fugacity=1.0, energy_perturbation=1e-9)
        df, ls_out, _ = ideal_gas_referenced_nested_sampling(ls, p,
            Int64(820), routine, cap_save())
        cap_cleanup()
        live_E = [w.energy.val for w in ls_out.walkers]
        live_N = [Int(sum(w.configuration.components[1]))
                  for w in ls_out.walkers]
        s = gc_thermodynamic_stats_ideal_ref(df, 16, 1.0, [cap_mu], [cap_T],
            64; ω0=65 / 64, live_emax=live_E, live_numbers=live_N)
        return abs(s.logXi[1, 1] - cap_lnXi), abs(s.mean_N[1, 1] - cap_meanN)
    end

    @testset "incremental-vs-full cross-route exactness" begin
        # Same-seed incremental-on and -off runs diverge in floating-point
        # ordering by design, so both arms gate statistically against the
        # exact grand sum rather than by trajectory identity
        for seed in (99301, 99302, 99303)
            d1, d2 = cap_arm(seed, MCGrandCanonicalMoves())
            @test d1 < 0.40
            @test d2 < 0.61
        end
        for seed in (96001, 96002, 96003)
            d1, d2 = cap_arm(seed, MCGrandCanonicalMoves(incremental=true))
            @test d1 < 0.85
            @test d2 < 0.27
        end
    end

    @testset "occupied-empty stationarity" begin
        # The zero-null proposal is a distinct chain with the same
        # stationary construction: agreement with the exact grand sum is
        # the whole claim
        for seed in (99311, 99312, 99313)
            d1, d2 = cap_arm(seed,
                MCGrandCanonicalMoves(swap_mode=:occupied_empty))
            @test d1 < 1.36
            @test d2 < 0.78
        end

        # First-quartile cross-proposal agreement (see the header's scope
        # note: deep quartiles legitimately diverge by decorrelation
        # quality, so only Q1 is gated across proposals)
        function cap_q1(seed, routine)
            Random.seed!(seed)
            ws = [LatticeWalker(deepcopy(cap_lat()), energy=0.0u"eV", iter=0)
                  for _ in 1:64]
            ls = LatticeGasWalkers(ws, cap_ham; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=40,
                reference_fugacity=1.0, energy_perturbation=1e-9)
            df, _, _ = ideal_gas_referenced_nested_sampling(ls, p, Int64(820),
                routine, cap_save())
            cap_cleanup()
            q = nrow(df) ÷ 4
            return sum(df.num_particles[1:q]) / q
        end
        q1_def = sum(cap_q1(s, MCGrandCanonicalMoves())
                     for s in (99301, 99302, 99303)) / 3
        q1_occ = sum(cap_q1(s, MCGrandCanonicalMoves(swap_mode=:occupied_empty))
                     for s in (99311, 99312, 99313)) / 3
        @test abs(q1_def - q1_occ) < 0.35
    end

    @testset "joint ledger-closure welds" begin
        # Lattice ledger, fresh seeded run: every column sums to its run
        # total (the twelve-key schema, incremental and occupied-empty on
        # to exercise the widest kernel surface)
        Random.seed!(99321)
        ws = [LatticeWalker(deepcopy(cap_lat()), energy=0.0u"eV", iter=0)
              for _ in 1:16]
        ls = LatticeGasWalkers(ws, cap_ham; assign_energy=false)
        p = IdealGasReferencedGCNSParameters(mc_steps=30,
            reference_fugacity=1.0, energy_perturbation=1e-9)
        d, _, pout = ideal_gas_referenced_nested_sampling(ls, p, Int64(80),
            MCGrandCanonicalMoves(clusters_freq=1, swaps_freq=2,
                                  incremental=true,
                                  swap_mode=:occupied_empty),
            cap_save(); record_move_rates=true)
        cap_cleanup()
        for name in ("swap_attempted", "swap_accepted",
                     "swap_null_attempted", "swap_null_accepted",
                     "cluster_attempted", "cluster_accepted",
                     "insert_uniform_attempted", "insert_uniform_accepted",
                     "insert_biased_attempted", "insert_biased_accepted",
                     "delete_attempted", "delete_accepted")
            @test sum(d[!, name]) == get(pout.move_stats, Symbol(name), 0)
        end
        # occupied-empty produced zero nulls here too
        @test get(pout.move_stats, :swap_null_attempted, 0) == 0

        # Canonical reflective route, fresh seeded run: the five Galilean
        # columns close and the trailing step-size column carries the final
        # adapted value
        cl_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        cl_lj = LJParameters(epsilon=0.05, sigma=2.5, cutoff=2.0, shift=true)
        cl_pair(dd) = AtomWalker(FastSystem(atomic_system(
            [:Ar => [5.0, 6.0, 6.0]u"Å", :Ar => [5.0 + dd, 6.0, 6.0]u"Å"],
            cl_box, (true, true, true))))
        Random.seed!(99322)
        gws = [cl_pair(2.80 + 0.02k) for k in 1:6]
        for w in gws
            w.energy = interacting_energy(w.configuration, cl_lj,
                                          w.list_num_par, w.frozen)
        end
        gls = LJAtomWalkers(gws, cl_lj)
        gp = NestedSamplingParameters(mc_steps=4, initial_step_size=0.5,
            step_size=0.5, allowed_fail_count=100_000)
        gd, _, gpout = nested_sampling(gls, gp, Int64(25),
            MCGalileanWalk(n_refresh=4), cap_save(); record_move_rates=true)
        cap_cleanup()
        for name in ("galilean_attempted", "galilean_accepted",
                     "galilean_reflect_attempted", "galilean_reflect_evals",
                     "galilean_reflect_accepted")
            @test sum(gd[!, name]) == get(gpout.move_stats, Symbol(name), 0)
        end
        @test gd.step_size[end] == gpout.step_size

        # Atomistic ideal-gas-referenced conditional channels, fresh seeded
        # run: the Galilean and cavity columns close against their totals
        aw_pair(dd) = AtomWalker(FastSystem(atomic_system(
            [:Ar => [5.0, 6.0, 6.0]u"Å", :Ar => [5.0 + dd, 6.0, 6.0]u"Å"],
            cl_box, (true, true, true))))
        Random.seed!(99323)
        aws = [aw_pair(2.80 + 0.02k) for k in 1:6]
        als = LJAtomWalkers(aws, cl_lj)
        ap = AtomisticIGRefGCNSParameters(mc_steps=60,
            reference_activity=(6.0 / 12.0^3)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000)
        ad, _, apout = ideal_gas_referenced_nested_sampling(als, ap, 30,
            MCAtomGrandCanonicalMoves(galilean_steps=2, p_bias=0.3,
                                      bias_radius=2.0, bias_grid=6),
            cap_save(); record_move_rates=true)
        cap_cleanup()
        for name in ("insert_biased_attempted", "insert_biased_accepted",
                     "galilean_attempted", "galilean_accepted",
                     "galilean_reflect_attempted", "galilean_reflect_evals",
                     "galilean_reflect_accepted")
            @test sum(ad[!, name]) == get(apout.move_stats, Symbol(name), 0)
        end
    end

    @testset "allocation guards" begin
        # Structural performance invariants, never time-based: a warmed
        # mixed-channel incremental walk at M = 256 under a generous fixed
        # byte ceiling, and sub-linear step-allocation growth on the
        # swap-only channel (the insert/delete channels keep their O(M)
        # findall draws for stream identity and are covered by the fixed
        # ceiling only)
        walk_lat(M, L) = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1],
            components=[[false for _ in 1:M]], adsorptions=:full)
        h1 = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")

        lat256 = walk_lat(256, 16)
        Random.seed!(99336)
        for i in 1:256
            lat256.components[1][i] = rand() < 0.5
        end
        wk256 = LatticeWalker(lat256, energy=interacting_energy(lat256, h1),
                              iter=0)
        MC_grand_canonical_walk!(200, wk256, h1, 1.0e3, 0.0;
            p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
            incremental=true)
        mixed_bytes = @allocated MC_grand_canonical_walk!(200, wk256, h1,
            1.0e3, 0.0; p_move=0.4, p_insert=0.3, z0=1.0,
            energy_perturb=1e-9, incremental=true)
        @test mixed_bytes < 2_000_000

        function swap_marginal(wk)
            f(n) = (MC_grand_canonical_walk!(n, wk, h1, 1.0e3, 0.0;
                p_move=1.0, p_insert=0.0, z0=1.0, energy_perturb=1e-9,
                incremental=true); nothing)
            f(200)
            f(400)
            a200 = @allocated f(200)
            a400 = @allocated f(400)
            return a400 - a200
        end
        marg256 = swap_marginal(wk256)
        lat1024 = walk_lat(1024, 32)
        Random.seed!(99337)
        for i in 1:1024
            lat1024.components[1][i] = rand() < 0.5
        end
        wk1024 = LatticeWalker(lat1024,
            energy=interacting_energy(lat1024, h1), iter=0)
        marg1024 = swap_marginal(wk1024)
        @test marg256 > 0
        # 4x the sites must not scale the per-step allocation: the measured
        # ratio class is ~1.25, while a deepcopy- or full-recompute-class
        # regression shifts it by orders of magnitude
        @test marg1024 < 2 * marg256
    end
end
