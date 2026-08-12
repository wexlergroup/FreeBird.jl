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
