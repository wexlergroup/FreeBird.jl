# Atomistic energy-sorted ideal-gas-referenced grand-canonical nested sampling.
#
# Calibration ledger (protocol identical to the shipped testsets):
# - Initializer moments at K = 200, z0V = 3, seeds {1,2,3,52521}: max devs
#   mean 0.25, var 0.339; gates ship at >= 3x (0.75 and 1.05).
# - Tie-block fixture (seed 52530): eviction charges are DETERMINISTIC
#   [log(5/6), log(4/5), log(3/4), log(2/3)] (exact logs of integer ratios,
#   architecture-exact). The refill outcome is trajectory-dependent and
#   architecture-sensitive at the ulp level (macOS-aarch64 refilled to
#   particle counts [2, 2, 2, 3, 3, 3]; ubuntu-x64 under Julia 1.12 to
#   [2, 2, 2, 3, 4, 5]: rounding drift flips accept/reject decisions inside
#   the refill walks), so the outcome asserts are trajectory-robust: the live
#   set is restored to size, every walker sits genuinely below the plateau
#   (negative energy, never rounding dust, thanks to the from-scratch
#   re-anchoring of accepted clones), and at least one refilled clone
#   re-entered at a different particle count (the N-changing refill contract).
# - End-to-end determinism at K = 16, seed 52522: two runs digit-identical.
@testset "atomistic ideal-gas-referenced GC-NS tests" begin
    using Random

    box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
    pbc = (true, true, true)
    seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
    mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                           empty(position(seed_at, :)), empty(species(seed_at, :)),
                           empty(mass(seed_at, :)))
    mkat(coords) = AtomWalker{1}(FastSystem(atomic_system([:Ar => SVector{3}(c)u"Å" for c in coords], box, pbc)))
    V = 1728.0
    lj0 = LJParameters(epsilon=0.0)
    routine = MCAtomGrandCanonicalMoves()
    mksave(tag) = SaveEveryN(df_filename="_igref_at_$(tag).csv",
                             wk_filename="_igref_at_$(tag).traj.extxyz",
                             ls_filename="_igref_at_$(tag).ls.extxyz",
                             n_traj=10^7, n_snap=10^7, n_info=10^7)
    clean(tag) = for f in ["_igref_at_$(tag).csv", "_igref_at_$(tag).traj.extxyz", "_igref_at_$(tag).ls.extxyz"]
        rm(f, force=true)
    end

    @testset "routine and parameter validation" begin
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(p_move=0.9, p_insert=0.2)
        @test_throws ArgumentError AtomisticIGRefGCNSParameters(reference_activity=0.0u"Å^-3")
        params = AtomisticIGRefGCNSParameters(reference_activity=(2.0 / V)u"Å^-3", species=:Ar)
        wf = AtomWalker{1}(mkempty())
        wf.frozen = [true]
        lsf = GenericAtomWalkers([wf], lj0; assign_energy=false)
        @test_throws ArgumentError FreeBird.SamplingSchemes._atomistic_igref_z0V(lsf, params)
        tilted = [[12.0, 0.0, 0.0], [2.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        wt = AtomWalker{1}(FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], tilted, pbc)))
        lst = GenericAtomWalkers([wt], lj0)
        @test_throws ArgumentError FreeBird.SamplingSchemes._atomistic_igref_z0V(lst, params)
        small = [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]u"Å"
        wm = AtomWalker{1}(FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], small, pbc)))
        lsm = GenericAtomWalkers([AtomWalker{1}(mkempty()), wm], lj0)
        @test_throws ArgumentError FreeBird.SamplingSchemes._atomistic_igref_z0V(lsm, params)
        # z0V folds the activity against the cell volume
        @test FreeBird.SamplingSchemes._atomistic_igref_z0V(
            GenericAtomWalkers([AtomWalker{1}(mkempty())], lj0), params) ≈ 2.0 atol=1e-12
    end

    @testset "initializer draws the reference law (seeded, calibrated)" begin
        Random.seed!(52521)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:200], lj0)
        for w in ls.walkers
            w.iter = 99
        end
        params = AtomisticIGRefGCNSParameters(reference_activity=(3.0 / V)u"Å^-3", species=:Ar)
        z0V = FreeBird.SamplingSchemes._atomistic_igref_z0V(ls, params)
        FreeBird.SamplingSchemes._init_atomistic_igref_walkers!(ls, params, z0V)
        ns = [w.list_num_par[1] for w in ls.walkers]
        @test abs(mean(ns) - 3.0) < 0.75
        @test abs(var(ns) - 3.0) < 1.05
        @test all(w.iter == 0 for w in ls.walkers)
        @test all(length(w.configuration) == w.list_num_par[1] for w in ls.walkers)
        @test all(w.energy == 0.0u"eV" for w in ls.walkers)
    end

    @testset "stall contract on an exactly degenerate liveset (deterministic)" begin
        Random.seed!(31312)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:8], lj0)
        params = AtomisticIGRefGCNSParameters(mc_steps=20, reference_activity=(4.0 / V)u"Å^-3",
                                              species=:Ar, allowed_fail_count=7)
        save = mksave("stall")
        df, ls_out, params_out = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ideal_gas_referenced_nested_sampling(ls, params, 500, routine, save)
        end
        clean("stall")
        @test nrow(df) == 0
        @test names(df) == ["iter", "emax", "num_particles", "log_compression"]
        @test length(ls_out.walkers) == 8
        @test all(w.energy == 0.0u"eV" for w in ls_out.walkers)
        @test params_out.fail_count == 7
        # warn-and-continue keeps the sibling contract: the loop runs its full budget
        Random.seed!(31313)
        ls2 = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:8], lj0)
        params2 = AtomisticIGRefGCNSParameters(mc_steps=20, reference_activity=(4.0 / V)u"Å^-3",
                                               species=:Ar, allowed_fail_count=7)
        save2 = mksave("stall2")
        df2, _, params2_out = @test_logs (:warn, r"Failed 7 times") match_mode=:any begin
            ideal_gas_referenced_nested_sampling(ls2, params2, 30, routine, save2; stop_on_stall=false)
        end
        clean("stall2")
        @test nrow(df2) == 0
        @test params2_out.fail_count < 7
    end

    @testset "plateau tie block: bit-exact charges and N-changing refill (seed 52530)" begin
        lj = LJParameters(epsilon=0.05, sigma=2.5, cutoff=2.5)
        rmin = 2.0^(1 / 6) * 2.5
        walkers = AtomWalker{1}[]
        push!(walkers, mkat([[1.0, 1.0, 1.0]]))
        push!(walkers, mkat([[1.0, 6.0, 1.0]]))
        push!(walkers, mkat([[6.0, 1.0, 1.0]]))
        push!(walkers, mkat([[6.0, 6.0, 6.0]]))
        push!(walkers, mkat([[3.0, 3.0, 3.0], [3.0 + rmin, 3.0, 3.0]]))
        push!(walkers, mkat([[9.0, 9.0, 9.0], [9.0 - rmin, 9.0, 9.0]]))
        ls = GenericAtomWalkers(walkers, lj)
        @test FreeBird.SamplingSchemes._tie_block_length(sort(ls.walkers, by=w -> w.energy, rev=true)) == 4
        Random.seed!(52530)
        params = AtomisticIGRefGCNSParameters(mc_steps=100, reference_activity=(2.0 / V)u"Å^-3",
                                              species=:Ar, allowed_fail_count=200, step_size=0.5)
        charges = Float64[]
        npars = Int[]
        for _ in 1:4
            iter, emax, n_par, ls, params, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(ls, params, routine)
            push!(charges, log_t)
            push!(npars, n_par)
        end
        @test charges == [log(5 / 6), log(4 / 5), log(3 / 4), log(2 / 3)]
        @test npars == [1, 1, 1, 1]
        @test length(ls.walkers) == 6
        @test params.plateau_refill_target == 0
        # trajectory-robust refill asserts (the exact particle-count multiset is
        # architecture-sensitive; see the calibration ledger)
        @test all(w.energy < 0.0u"eV" for w in ls.walkers)
        @test any(w.list_num_par[1] != 2 for w in ls.walkers)
        @test all(1 <= w.list_num_par[1] <= 8 for w in ls.walkers)
        @test all(length(w.configuration) == w.list_num_par[1] for w in ls.walkers)
    end

    @testset "dilute interacting end-to-end with same-seed determinism" begin
        function run_e2e(seed)
            Random.seed!(seed)
            ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:16],
                                    LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5))
            params = AtomisticIGRefGCNSParameters(mc_steps=60, reference_activity=(6.0 / V)u"Å^-3",
                                                  species=:Ar, allowed_fail_count=50)
            save = mksave("e2e")
            df, lso, pout = ideal_gas_referenced_nested_sampling(ls, params, 120, routine, save)
            clean("e2e")
            return df, [ustrip(u"eV", w.energy) for w in lso.walkers], pout
        end
        df1, live1, pout1 = run_e2e(52522)
        df2, live2, _ = run_e2e(52522)
        @test df1 == df2
        @test live1 == live2
        @test nrow(df1) == 120
        @test names(df1) == ["iter", "emax", "num_particles", "log_compression"]
        @test issorted(df1.emax, rev=true)
        @test all(df1.log_compression .< 0.0)
        @test all(df1.num_particles .>= 0)
        @test haskey(pout1.move_stats, :insert_attempted)
        @test haskey(pout1.move_stats, :delete_attempted)
        @test pout1.plateau_refill_target == 0
    end

    @testset "parameter reuse: dirty runtime state cannot leak into a run" begin
        function run_with(params_init)
            Random.seed!(52523)
            ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:8],
                                    LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5))
            save = mksave("reuse")
            df, lso, pout = ideal_gas_referenced_nested_sampling(ls, params_init, 60, routine, save)
            clean("reuse")
            return df, [ustrip(u"eV", w.energy) for w in lso.walkers], pout
        end
        clean_params() = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / V)u"Å^-3",
                                                      species=:Ar, allowed_fail_count=50)
        dirty = clean_params()
        dirty.plateau_refill_target = 3
        dirty.fail_count = 5
        dirty.move_stats[:stale] = 99
        df_c, live_c, _ = run_with(clean_params())
        df_d, live_d, pout_d = run_with(dirty)
        @test df_c == df_d
        @test live_c == live_d
        @test !haskey(pout_d.move_stats, :stale)
    end

    @testset "default-cadence saves skip zero-atom frames instead of crashing" begin
        Random.seed!(52541)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:8], LJParameters(epsilon=0.0))
        params = AtomisticIGRefGCNSParameters(mc_steps=10, reference_activity=(0.5 / V)u"Å^-3",
                                              species=:Ar, allowed_fail_count=5)
        save = SaveEveryN(df_filename="_igref_at_guard.csv",
                          wk_filename="_igref_at_guard.traj.extxyz",
                          ls_filename="_igref_at_guard.ls.extxyz",
                          n_traj=1, n_snap=1, n_info=10^7)
        df, lso, _ = @test_logs (:warn, r"zero-atom frame") match_mode=:any begin
            ideal_gas_referenced_nested_sampling(ls, params, 20, routine, save)
        end
        for f in ["_igref_at_guard.csv", "_igref_at_guard.traj.extxyz", "_igref_at_guard.ls.extxyz"]
            rm(f, force=true)
        end
        @test length(lso.walkers) == 8
    end

    @testset "observables and dead-point callback handle variable N" begin
        Random.seed!(52540)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:8],
                                LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5))
        params = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / V)u"Å^-3",
                                              species=:Ar, allowed_fail_count=50)
        save = mksave("obs")
        seen_ns = Int[]
        df, _, _ = ideal_gas_referenced_nested_sampling(
            ls, params, 40, routine, save;
            observables=[:n_obs => cfg -> Float64(length(cfg))],
            dead_point_callback=(iter, walker) -> push!(seen_ns, walker.list_num_par[1]))
        clean("obs")
        @test "n_obs" in names(df)
        @test df.n_obs == Float64.(df.num_particles)
        @test seen_ns == df.num_particles
    end
end

@testset "Interacting-regime controls" begin
    using Random
    # Fixtures calibrated in-session; see the assertions for the measured
    # separations. Bound-dimer walkers at slightly different separations give a
    # continuous energy ladder with collapsed insertion/deletion acceptance at
    # tiny z0V (Metropolis ratio z0V/(N+1)), while displacements stay healthy.
    hard_L = 12.0
    hard_box = [[hard_L * u"Å", 0u"Å", 0u"Å"],
                [0u"Å", hard_L * u"Å", 0u"Å"],
                [0u"Å", 0u"Å", hard_L * u"Å"]]
    hard_lj = LJParameters(epsilon=0.05, sigma=2.5, cutoff=2.0, shift=true)
    hard_pair(r) = AtomWalker(FastSystem(periodic_system(
        [:Ar => [0.3, 0.5, 0.5], :Ar => [0.3 + r / hard_L, 0.5, 0.5]],
        hard_box, fractional=true)))

    @testset "Routine and parameter validation" begin
        r = MCAtomGrandCanonicalMoves()
        @test r.step_rate_source == :mixed
        @test r.mc_steps_per_particle == 0.0
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(step_rate_source=:invalid)
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(mc_steps_per_particle=-1.0)
        p = AtomisticIGRefGCNSParameters()
        @test p.allowed_fail_count == 100
        @test p.refill_fail_budget == 0
        @test_throws ArgumentError AtomisticIGRefGCNSParameters(refill_fail_budget=-1)
    end

    @testset "Walk-length arithmetic" begin
        p = AtomisticIGRefGCNSParameters(mc_steps=100)
        @test FreeBird.SamplingSchemes._gc_walk_length(p, MCAtomGrandCanonicalMoves(), 50) == 100
        @test FreeBird.SamplingSchemes._gc_walk_length(p, MCAtomGrandCanonicalMoves(mc_steps_per_particle=2.0), 50) == 200
        @test FreeBird.SamplingSchemes._gc_walk_length(p, MCAtomGrandCanonicalMoves(mc_steps_per_particle=1.5), 2) == 103
    end

    @testset "Step-rate source A/B" begin
        # Calibrated at seed 3000: the mixed rate ~ 0.3*A_move + 0.35*A_ins +
        # 0.35*A_del sits below the 0.2 band edge (measured A_ins ~ A_del ~ 0.010,
        # A_move 0.56-0.64),
        # so :mixed shrinks the step (measured 0.028) while :move holds it
        # inside the wide band (measured 0.104)
        results = Dict{Symbol,Float64}()
        stats = Dict{Symbol,Dict{Symbol,Int}}()
        for src in (:mixed, :move)
            Random.seed!(3000)
            ls = LJAtomWalkers([hard_pair(2.80 + 0.02k) for k in 1:6], hard_lj)
            params = AtomisticIGRefGCNSParameters(mc_steps=100,
                reference_activity=(0.05 / hard_L^3)u"Å^-3", species=:Ar,
                step_size=0.1, accept_range=(0.2, 0.9), allowed_fail_count=100_000)
            routine = MCAtomGrandCanonicalMoves(p_move=0.3, p_insert=0.35, step_rate_source=src)
            z0V = FreeBird.SamplingSchemes._atomistic_igref_z0V(ls, params)
            for i in 1:30
                FreeBird.SamplingSchemes.nested_sampling_step!(ls, params, routine; ns_iteration=i, z0V=z0V)
            end
            results[src] = params.step_size
            stats[src] = copy(params.move_stats)
        end
        for src in (:mixed, :move)
            ms = stats[src]
            @test ms[:move_accepted] / ms[:move_attempted] > 0.5
            @test ms[:insert_accepted] / ms[:insert_attempted] < 0.05
        end
        @test results[:mixed] < 0.05
        @test results[:move] > 0.09
    end

    @testset "Refill failure budget" begin
        # Bit-exact duplicate dimers form a tie block above three slightly
        # deeper dimers; five-step walks at step 1.5 rarely land below the tie
        # energy, so the cumulative budget (allowed_fail_count = 2) abandons
        # the refill where a dedicated budget completes it. Calibrated at
        # seeds 4000/4001: n_live 3 vs 6 at both.
        for (budget, expected) in ((0, 3), (300, 6))
            Random.seed!(4000)
            dup = hard_pair(2.90)
            ls = LJAtomWalkers([deepcopy(dup), deepcopy(dup), deepcopy(dup),
                                hard_pair(2.82), hard_pair(2.83), hard_pair(2.84)], hard_lj)
            params = AtomisticIGRefGCNSParameters(mc_steps=5,
                reference_activity=(0.05 / hard_L^3)u"Å^-3", species=:Ar,
                step_size=1.5, allowed_fail_count=2, refill_fail_budget=budget)
            routine = MCAtomGrandCanonicalMoves(p_move=0.6, p_insert=0.2)
            z0V = FreeBird.SamplingSchemes._atomistic_igref_z0V(ls, params)
            for i in 1:6
                FreeBird.SamplingSchemes.nested_sampling_step!(ls, params, routine; ns_iteration=i, z0V=z0V)
            end
            @test length(ls.walkers) == expected
        end
    end

    @testset "Stall-continue resets the step size" begin
        # Zero-interaction degenerate live set: every replacement fails, and on
        # the stop_on_stall=false branch the step size returns to initial
        Random.seed!(5000)
        seeds = [AtomWalker(FastSystem(periodic_system(
            [:Ar => [0.5, 0.5, 0.5]], hard_box, fractional=true))) for _ in 1:4]
        ls = GenericAtomWalkers(seeds, IdealGasParameters())
        params = AtomisticIGRefGCNSParameters(mc_steps=10,
            reference_activity=(2.0 / hard_L^3)u"Å^-3", species=:Ar,
            initial_step_size=0.5, step_size=1.7, allowed_fail_count=3)
        save = SaveEveryN(df_filename="_t_stall.csv", wk_filename="_t_stall.traj.extxyz",
                          ls_filename="_t_stall.ls.extxyz",
                          n_traj=10^7, n_snap=10^7, n_info=10^7)
        df, ls, params = ideal_gas_referenced_nested_sampling(
            ls, params, 3, MCAtomGrandCanonicalMoves(), save; stop_on_stall=false)
        for f in ("_t_stall.csv", "_t_stall.traj.extxyz", "_t_stall.ls.extxyz")
            rm(f, force=true)
        end
        @test params.step_size == params.initial_step_size
        @test nrow(df) == 0
    end

    @testset "Per-iteration acceptance ledger" begin
        Random.seed!(6000)
        ls = LJAtomWalkers([hard_pair(2.80 + 0.02k) for k in 1:6], hard_lj)
        params = AtomisticIGRefGCNSParameters(mc_steps=60,
            reference_activity=(6.0 / hard_L^3)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000)
        save = SaveEveryN(df_filename="_t_rates.csv", wk_filename="_t_rates.traj.extxyz",
                          ls_filename="_t_rates.ls.extxyz",
                          n_traj=10^7, n_snap=10^7, n_info=10^7)
        df, ls, params = ideal_gas_referenced_nested_sampling(
            ls, params, 40, MCAtomGrandCanonicalMoves(), save; record_move_rates=true)
        rate_names = ["move_attempted", "move_accepted", "insert_attempted",
                      "insert_accepted", "delete_attempted", "delete_accepted"]
        @test names(df) == vcat(["iter", "emax", "num_particles", "log_compression"], rate_names)
        # Closure: recorded per-iteration deltas sum to the run totals
        for name in rate_names
            @test sum(df[!, name]) == get(params.move_stats, Symbol(name), 0)
        end
        # The default keeps the shipped four-column schema
        Random.seed!(6000)
        ls2 = LJAtomWalkers([hard_pair(2.80 + 0.02k) for k in 1:6], hard_lj)
        params2 = AtomisticIGRefGCNSParameters(mc_steps=60,
            reference_activity=(6.0 / hard_L^3)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000)
        save2 = SaveEveryN(df_filename="_t_rates.csv", wk_filename="_t_rates.traj.extxyz",
                           ls_filename="_t_rates.ls.extxyz",
                           n_traj=10^7, n_snap=10^7, n_info=10^7)
        df2, _, _ = ideal_gas_referenced_nested_sampling(
            ls2, params2, 40, MCAtomGrandCanonicalMoves(), save2)
        @test names(df2) == ["iter", "emax", "num_particles", "log_compression"]
        # Same seed, recording on or off: identical sampling trajectory
        @test df2.emax == df.emax
        for f in ("_t_rates.csv", "_t_rates.traj.extxyz", "_t_rates.ls.extxyz")
            rm(f, force=true)
        end
        # The rate columns are reserved: observables cannot collide with them
        @test_throws ArgumentError FreeBird.SamplingSchemes._validate_observables(
            [:insert_accepted => cfg -> 0.0], ls)
    end

    @testset "Minimum-image cutoff guard" begin
        cfg = hard_pair(2.9).configuration
        over = LJParameters(epsilon=0.05, sigma=2.5, cutoff=2.5, shift=true)   # 6.25 > 6.0
        under = LJParameters(epsilon=0.05, sigma=2.5, cutoff=2.0, shift=true)  # 5.0 <= 6.0
        untr = LJParameters(epsilon=0.05, sigma=2.5)                           # Inf: never warns
        @test_logs (:warn, r"minimum-image") FreeBird.SamplingSchemes._warn_min_image_cutoff(over, cfg)
        @test_logs FreeBird.SamplingSchemes._warn_min_image_cutoff(under, cfg)
        @test_logs FreeBird.SamplingSchemes._warn_min_image_cutoff(untr, cfg)
        @test SamplingSchemes.AbstractPotentials._max_interaction_range(IdealGasParameters()) === missing
    end
end

@testset "Galilean routine and grand-canonical alternation" begin
    using Random
    # Calibration ledger (gates at >= 3x the max three-seed deviation, per stat):
    # - Canonical cross-route log_Z_N (N = 3 LJ-fluid sector, 12 A box, K = 32,
    #   500 steps, reduced at 200 K; RW seeds 7100x vs Galilean seeds 7110x):
    #   diffs 0.0901, 0.0776, 0.0349; gate 0.28 nats.
    # - Ideal-gas-referenced alternation on/off (dilute LJ, z0V = 6, K = 16,
    #   150 steps, logXi at the 300 K anchor; seeds 7200x/7210x): diffs 0.0322,
    #   0.0184, 0.1709; gate 0.55 nats.
    gal_L = 12.0
    gal_box = [[gal_L * u"Å", 0u"Å", 0u"Å"],
               [0u"Å", gal_L * u"Å", 0u"Å"],
               [0u"Å", 0u"Å", gal_L * u"Å"]]
    gal_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=3.0, shift=true)
    function gal_uwalker(N)
        while true
            coor = [:Ar => [rand(), rand(), rand()] for _ in 1:N]
            sys = FastSystem(periodic_system(coor, gal_box, fractional=true))
            E = ustrip(u"eV", interacting_energy(sys, gal_lj))
            (isfinite(E) && E < 100.0) && return AtomWalker(sys)
        end
    end
    gal_save(tag) = SaveEveryN(df_filename="_$(tag).csv", wk_filename="_$(tag).w.extxyz",
                               ls_filename="_$(tag).l.extxyz",
                               n_traj=10^8, n_snap=10^8, n_info=10^8)
    gal_rm(tag) = for f in ("_$(tag).csv", "_$(tag).w.extxyz", "_$(tag).l.extxyz")
        rm(f, force=true)
    end

    @testset "routine and field validation" begin
        @test MCGalileanWalk().n_refresh == 8
        @test_throws ArgumentError MCGalileanWalk(n_refresh=0)
        r = MCAtomGrandCanonicalMoves()
        @test r.galilean_steps == 0 && r.galilean_n_refresh == 8 && r.galilean_step_size == 0.5
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(galilean_steps=-1)
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(galilean_step_size=0.0)
    end

    @testset "canonical cross-route agreement (seeded, calibrated)" begin
        function run_route(routine, seed, tag)
            Random.seed!(seed)
            ls = LJAtomWalkers([gal_uwalker(3) for _ in 1:32], gal_lj)
            p = NestedSamplingParameters(mc_steps=300, initial_step_size=0.3, step_size=0.3,
                step_size_lo=0.01, step_size_up=2.0, accept_range=(0.25, 0.75),
                allowed_fail_count=1000, energy_perturbation=1e-12)
            df, fls, _ = nested_sampling(ls, p, 500, routine, gal_save(tag))
            gal_rm(tag)
            live = [ustrip(u"eV", w.energy) for w in fls.walkers]
            out = gc_thermodynamic_stats_fixed_N(
                [DataFrame(iter=Int[], emax=Float64[]), df], [0, 3],
                (gal_L^3)u"Å^3", 39.948u"u", [-0.20u"eV"], [200.0]u"K";
                n_walkers=32, live_emax=[Float64[], live])
            return out.log_Z_N[2, 1], df
        end
        a, df_rw = run_route(MCRandomWalkClone(), 71001, "xrw")
        b, df_gw = run_route(MCGalileanWalk(n_refresh=6), 71101, "xgw")
        @test abs(a - b) < 0.28
        # the Galilean route drives a genuine descent with the serial ledger schema
        @test names(df_gw) == ["iter", "emax", "log_compression"]
        @test nrow(df_gw) >= 400
        @test issorted(df_gw.emax, rev=true) || maximum(diff(df_gw.emax)) <= 1e-12
    end

    @testset "step-method contract and plateau cooperation (tie fixture)" begin
        # bit-exact duplicate dimers above two deeper dimers: the new step
        # method's plateau branch must charge the deterministic Fowlie-Handley-Su
        # schedule (exact logs of integer ratios) and refill through the
        # reflective kernel
        gal_pair(r) = AtomWalker(FastSystem(periodic_system(
            [:Ar => [0.3, 0.5, 0.5], :Ar => [0.3 + r / gal_L, 0.5, 0.5]],
            gal_box, fractional=true)))
        Random.seed!(76543)
        dup = gal_pair(2.90)
        ls = LJAtomWalkers([deepcopy(dup), deepcopy(dup), deepcopy(dup),
                            gal_pair(2.82), gal_pair(2.84)], gal_lj)
        p = NestedSamplingParameters(mc_steps=5, initial_step_size=0.5, step_size=0.5,
            step_size_lo=0.01, step_size_up=2.0, accept_range=(0.25, 0.75),
            allowed_fail_count=1000, energy_perturbation=1e-12)
        df, fls, _ = nested_sampling(ls, p, 4, MCGalileanWalk(n_refresh=4), gal_save("tie"))
        gal_rm("tie")
        @test df.log_compression[1:3] ≈ [log(4 / 5), log(3 / 4), log(2 / 3)] atol = 1e-14
        @test length(fls.walkers) >= 2
        # a direct step call returns the serial steps' five-value shape
        Random.seed!(76544)
        ls2 = LJAtomWalkers([gal_pair(2.80 + 0.02k) for k in 1:4], gal_lj)
        p2 = NestedSamplingParameters(mc_steps=5, initial_step_size=0.5, step_size=0.5,
            step_size_lo=0.01, step_size_up=2.0, accept_range=(0.25, 0.75),
            allowed_fail_count=1000, energy_perturbation=1e-12)
        out = FreeBird.SamplingSchemes.nested_sampling_step!(ls2, p2, MCGalileanWalk())
        @test length(out) == 5
        @test out[3] === ls2
        @test out[5] isa Union{Missing,Float64}
    end

    @testset "grand-canonical alternation (seeded, calibrated)" begin
        kb = 8.617333262e-5
        lam300 = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(39.948u"u", 300.0u"K"))
        mu300 = (kb * 300 * (log(6.0 / gal_L^3) + 3 * log(lam300))) * u"eV"
        function run_igref(routine, seed, tag)
            Random.seed!(seed)
            ls = LJAtomWalkers([gal_uwalker(1) for _ in 1:16], gal_lj)
            params = AtomisticIGRefGCNSParameters(mc_steps=100,
                reference_activity=(6.0 / gal_L^3)u"Å^-3", species=:Ar,
                allowed_fail_count=100_000)
            df, fls, _ = ideal_gas_referenced_nested_sampling(ls, params, 150, routine, gal_save(tag))
            gal_rm(tag)
            live_e = [ustrip(u"eV", w.energy) for w in fls.walkers]
            live_n = [w.list_num_par[1] for w in fls.walkers]
            out = gc_thermodynamic_stats_ideal_ref(df, (gal_L^3)u"Å^3", 39.948u"u",
                (6.0 / gal_L^3)u"Å^-3", [mu300], [300.0]u"K";
                live_emax=live_e, live_numbers=live_n)
            return out.logXi[1, 1], df
        end
        off, df_off = run_igref(MCAtomGrandCanonicalMoves(), 72001, "aoff")
        on, df_on = run_igref(MCAtomGrandCanonicalMoves(galilean_steps=3, galilean_step_size=1.0),
                              72101, "aon")
        @test abs(off - on) < 0.55
        # the channel changes nothing structural: same schema and stall contract
        @test names(df_on) == names(df_off)
        # zero-default draw-count identity: an explicit galilean_steps = 0 run
        # reproduces the default routine digit for digit
        _, e1 = run_igref(MCAtomGrandCanonicalMoves(), 73001, "ab1")
        _, e2 = run_igref(MCAtomGrandCanonicalMoves(galilean_steps=0), 73001, "ab2")
        @test e1.emax == e2.emax
        @test e1.num_particles == e2.num_particles
    end
end

@testset "cavity-bias plumbing for the ideal-gas-referenced driver" begin
    using Random
    # Calibration ledger: driver-level biased (p_bias = 0.4, r = 2.3 A,
    # grid = 10) vs uniform runs on the dilute LJ fixture (z0V = 6, K = 16,
    # 150 steps, logXi at the 300 K anchor; seeds 7400x/7410x): diffs 0.2027,
    # 0.6164, 0.1150, mixed signs; gate 1.85 nats (3x the max). This wiring
    # gate is deliberately loose: the kernel-level stationarity, invariance,
    # and locator-oracle tests carry the sharp end of the bias correctness.
    cav_L = 12.0
    cav_box = [[cav_L * u"Å", 0u"Å", 0u"Å"],
               [0u"Å", cav_L * u"Å", 0u"Å"],
               [0u"Å", 0u"Å", cav_L * u"Å"]]
    cav_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=3.0, shift=true)
    function cav_uwalker(N)
        while true
            coor = [:Ar => [rand(), rand(), rand()] for _ in 1:N]
            sys = FastSystem(periodic_system(coor, cav_box, fractional=true))
            E = ustrip(u"eV", interacting_energy(sys, cav_lj))
            (isfinite(E) && E < 100.0) && return AtomWalker(sys)
        end
    end
    cav_save(tag) = SaveEveryN(df_filename="_$(tag).csv", wk_filename="_$(tag).w.extxyz",
                               ls_filename="_$(tag).l.extxyz",
                               n_traj=10^8, n_snap=10^8, n_info=10^8)
    cav_rm(tag) = for f in ("_$(tag).csv", "_$(tag).w.extxyz", "_$(tag).l.extxyz")
        rm(f, force=true)
    end
    function cav_run(routine, seed, tag; n_steps=150)
        Random.seed!(seed)
        ls = LJAtomWalkers([cav_uwalker(1) for _ in 1:16], cav_lj)
        params = AtomisticIGRefGCNSParameters(mc_steps=100,
            reference_activity=(6.0 / cav_L^3)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000)
        df, fls, pout = ideal_gas_referenced_nested_sampling(ls, params, n_steps, routine, cav_save(tag))
        cav_rm(tag)
        return df, fls, pout
    end

    @testset "routine validation" begin
        r = MCAtomGrandCanonicalMoves()
        @test r.p_bias == 0.0 && r.bias_radius == 0.0 && r.bias_grid == 0
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(p_bias=1.5, bias_radius=2.0, bias_grid=8)
        @test_throws ArgumentError MCAtomGrandCanonicalMoves(p_bias=0.5)
        @test_logs (:warn, r"p_bias = 1") match_mode = :any begin
            MCAtomGrandCanonicalMoves(p_bias=1.0, bias_radius=2.0, bias_grid=8)
        end
    end

    @testset "defaults stream identity" begin
        df1, _, _ = cav_run(MCAtomGrandCanonicalMoves(), 75001, "cd1")
        df2, _, _ = cav_run(MCAtomGrandCanonicalMoves(p_bias=0.0, bias_radius=0.0, bias_grid=0),
                            75001, "cd2")
        @test df1.emax == df2.emax
        @test df1.num_particles == df2.num_particles
    end

    @testset "driver-level biased run (seeded, calibrated)" begin
        kb = 8.617333262e-5
        lam300 = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(39.948u"u", 300.0u"K"))
        mu300 = (kb * 300 * (log(6.0 / cav_L^3) + 3 * log(lam300))) * u"eV"
        function reduce_run(df, fls)
            live_e = [ustrip(u"eV", w.energy) for w in fls.walkers]
            live_n = [w.list_num_par[1] for w in fls.walkers]
            out = gc_thermodynamic_stats_ideal_ref(df, (cav_L^3)u"Å^3", 39.948u"u",
                (6.0 / cav_L^3)u"Å^-3", [mu300], [300.0]u"K";
                live_emax=live_e, live_numbers=live_n)
            return out.logXi[1, 1]
        end
        df_off, fls_off, _ = cav_run(MCAtomGrandCanonicalMoves(), 74001, "coff")
        df_on, fls_on, p_on = cav_run(MCAtomGrandCanonicalMoves(p_bias=0.4, bias_radius=2.3, bias_grid=10),
                                      74101, "con")
        @test names(df_on) == names(df_off)
        @test abs(reduce_run(df_off, fls_off) - reduce_run(df_on, fls_on)) < 1.85
        # the biased sub-channel fired and its counters flowed into the run totals
        @test get(p_on.move_stats, :insert_biased_attempted, 0) > 0
        @test get(p_on.move_stats, :insert_biased_accepted, 0) <=
              get(p_on.move_stats, :insert_biased_attempted, 0)
    end

    @testset "refill-path forwarding" begin
        # the shipped tie fixture with the channel enabled: the deterministic
        # eviction charges are unchanged (exact logs of integer ratios), and the
        # refill path runs the biased kernel without disturbing the contract
        Random.seed!(52530)
        dup_sys = FastSystem(periodic_system([:Ar => [0.30, 0.5, 0.5]], cav_box, fractional=true))
        dups = [AtomWalker{1}(deepcopy(dup_sys)) for _ in 1:4]
        for d in dups
            d.energy = 0.0u"eV"
        end
        deeper = [cav_uwalker(2) for _ in 1:2]
        ls = LJAtomWalkers(vcat(dups, deeper), cav_lj)
        params = AtomisticIGRefGCNSParameters(mc_steps=40,
            reference_activity=(2.0 / cav_L^3)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000)
        routine = MCAtomGrandCanonicalMoves(p_bias=0.4, bias_radius=2.3, bias_grid=10)
        df, fls, _ = ideal_gas_referenced_nested_sampling(ls, params, 30, routine, cav_save("crf"))
        cav_rm("crf")
        # walkers below zero exist (the deeper pair walkers), so the E = 0
        # plateau of the four duplicates is a genuine tie block
        charges = df.log_compression[1:3]
        @test charges ≈ [log(5 / 6), log(4 / 5), log(3 / 4)] atol = 1e-14
        @test length(fls.walkers) >= 2
    end
end
