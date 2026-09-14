# Chemical-potential ordering of the atomistic ideal-gas-referenced route
# (the `chemical_potential` field of `AtomisticIGRefGCNSParameters`; the step
# sorts, detects plateaus, selects parents, and accepts on Ω = E − μN, the
# ledger gains an :omega column, and the reduction is unchanged).
#
# Calibration ledger (protocol identical to the shipped testsets: gates ship at
# >= 3x the maximum three-seed deviation, stated per statistic):
# - Kernel stationarity under the Ω ceiling (150 independent walkers x 400
#   steps, ideal gas, z0V = 3, mu = -1 eV, ceiling 3.5 eV, so the stationary
#   law is Poisson(3) conditioned on N <= 3 with mean 51/26 and variance
#   0.88314; seeds {1,2,3,86001}): max devs mean 0.128 (at the shipped seed),
#   var 0.137 (seed 3); gates ship at >= 3x (0.39 and 0.42).
#   Under mu = +1 eV, ceiling -2.5 eV, n_max = 6 (Poisson(3) conditioned on
#   3 <= N <= 6, mean 3.9588, variance 0.9674; seeds {1,2,3,86002}): max devs
#   mean 0.175, var 0.196 (both at the shipped seed); gates ship at >= 3x
#   (0.53 and 0.59).
# - Ideal-gas descent under Ω ordering (K = 256, z0V = 4, mu = -0.02 eV,
#   unbounded; the reduction at the reference activity is exact, and the
#   activity-linearity, occupancy-moment, and p_N gates are calibrated over
#   seeds {1,2,3,86003}): max devs logXi 0.077 (seed 2), mean_N 0.413 (seed 2),
#   var_N 1.58 and total variation 0.0933 (both at the shipped seed; the
#   total-variation gate clears 3x by 2e-5); gates ship at >= 3x (0.24, 1.25,
#   4.8, 0.28) and also bound the bounded (n_max = 6) and
#   the mu > 0 (n_max = 8) descents, whose shipped-seed deviations are below
#   0.031.
# - Every other testset is an exact contract (bit-identity, closed forms,
#   ordering invariants, schema) and needs no calibration.
@testset "chemical-potential ordering of the atomistic ideal-gas-referenced route" begin
    using Random

    om_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
    om_pbc = (true, true, true)
    om_seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], om_box, om_pbc))
    om_mkempty() = FastSystem(cell_vectors(om_seed_at), periodicity(om_seed_at),
                              empty(position(om_seed_at, :)), empty(species(om_seed_at, :)),
                              empty(mass(om_seed_at, :)))
    om_V = 1728.0
    om_Vq = 1728.0u"Å^3"
    # cutoff 2.0 sigma = 5 A stays inside the half cell, so the driver's
    # minimum-image warning stays silent
    om_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.0)
    om_lj0 = LJParameters(epsilon=0.0)
    om_kb = 8.617333262e-5
    om_mass = 39.948u"u"
    om_lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(om_mass, T * u"K"))
    om_mu_for(zV, T) = (om_kb * T * (log(zV / om_V) + 3 * log(om_lam(T)))) * u"eV"
    om_logfact(n) = n == 0 ? 0.0 : sum(log, 1:n)
    om_lpps = FreeBird.AnalysisTools._log_poisson_partial_sum
    om_save(tag) = SaveEveryN(df_filename="_igomega_$(tag).csv",
                              wk_filename="_igomega_$(tag).traj.extxyz",
                              ls_filename="_igomega_$(tag).ls.extxyz",
                              n_traj=10^7, n_snap=10^7, n_info=10^7)
    om_clean(tag) = for f in ["_igomega_$(tag).csv", "_igomega_$(tag).traj.extxyz",
                              "_igomega_$(tag).ls.extxyz"]
        rm(f, force=true)
    end
    function om_walker(n; species=:Ar, box=(12.0, 12.0, 12.0), mk=om_mkempty)
        w = AtomWalker{1}(mk())
        for _ in 1:n
            pos = SVector(rand() * box[1], rand() * box[2], rand() * box[3])u"Å"
            FreeBird.AbstractWalkers.insert_particle!(w, pos, species)
        end
        return w
    end

    # surface fixture of the surface-route testsets
    om_sbox = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
    om_spbc = (true, true, false)
    om_sV = 1500.0
    om_surf_sys = FastSystem(atomic_system(
        [:H => [2.5, 2.5, 2.0]u"Å", :H => [7.5, 2.5, 2.0]u"Å",
         :H => [2.5, 7.5, 2.0]u"Å", :H => [7.5, 7.5, 2.0]u"Å"], om_sbox, om_spbc))
    om_mksurf() = AtomWalker(deepcopy(om_surf_sys); freeze_species=[:H])
    om_smkempty() = FastSystem(cell_vectors(om_surf_sys), periodicity(om_surf_sys),
                               empty(position(om_surf_sys, :)), empty(species(om_surf_sys, :)),
                               empty(mass(om_surf_sys, :)))
    om_cps = CompositeParameterSets(2, [LJParameters(epsilon=0.001, sigma=2.5, cutoff=1.8, shift=true),
                                        LJParameters(epsilon=0.003, sigma=2.5, cutoff=1.8, shift=true),
                                        LJParameters(epsilon=0.01, sigma=2.5, cutoff=1.8, shift=true)])
    om_cps0 = CompositeParameterSets(2, [om_lj0, om_lj0, om_lj0])
    om_sliveset(K) = LJSurfaceWalkers([AtomWalker{1}(om_smkempty()) for _ in 1:K], om_cps, om_mksurf();
                                      assign_energy=true)

    @testset "constructor, field, and ledger schema" begin
        @test AtomisticIGRefGCNSParameters().chemical_potential == 0.0u"eV"
        @test AtomisticIGRefGCNSParameters(chemical_potential=-0.02u"eV").chemical_potential == -0.02u"eV"
        @test AtomisticIGRefGCNSParameters(chemical_potential=0.5u"eV").chemical_potential == 0.5u"eV"
        # a unitless chemical potential is a loud type error, never a silent eV
        @test_throws TypeError AtomisticIGRefGCNSParameters(chemical_potential=0.01)
        # a NaN would make every ceiling comparison false; a negative zero is
        # normalized to the default's exact +0.0
        @test_throws ArgumentError AtomisticIGRefGCNSParameters(chemical_potential=NaN * u"eV")
        @test !signbit(ustrip(u"eV", AtomisticIGRefGCNSParameters(chemical_potential=-0.0u"eV").chemical_potential))
        # the :omega column is reserved
        Random.seed!(86010)
        ls = GenericAtomWalkers([AtomWalker{1}(om_mkempty()) for _ in 1:4], om_lj)
        p = AtomisticIGRefGCNSParameters(mc_steps=10, reference_activity=(2.0 / om_V)u"Å^-3",
                                         species=:Ar, allowed_fail_count=20,
                                         chemical_potential=-0.01u"eV")
        @test_throws ArgumentError ideal_gas_referenced_nested_sampling(
            ls, p, 2, MCAtomGrandCanonicalMoves(), om_save("res");
            observables=[:omega => cfg -> 0.0])
        om_clean("res")
        # schema: default unchanged, nonzero mu appends :omega after the compression
        # column and before the acceptance columns and observables
        function schema_run(mu; record=false, obs=nothing)
            Random.seed!(86011)
            lsr = GenericAtomWalkers([AtomWalker{1}(om_mkempty()) for _ in 1:8], om_lj)
            pr = AtomisticIGRefGCNSParameters(mc_steps=20, reference_activity=(3.0 / om_V)u"Å^-3",
                                              species=:Ar, allowed_fail_count=50,
                                              chemical_potential=mu)
            df, _, _ = ideal_gas_referenced_nested_sampling(
                lsr, pr, 12, MCAtomGrandCanonicalMoves(), om_save("sch");
                record_move_rates=record, observables=obs)
            om_clean("sch")
            return df
        end
        @test names(schema_run(0.0u"eV")) == ["iter", "emax", "num_particles", "log_compression"]
        @test names(schema_run(-0.01u"eV")) == ["iter", "emax", "num_particles", "log_compression", "omega"]
        @test names(schema_run(-0.01u"eV"; record=true)) ==
              ["iter", "emax", "num_particles", "log_compression", "omega",
               "move_attempted", "move_accepted", "insert_attempted", "insert_accepted",
               "delete_attempted", "delete_accepted", "step_size"]
        @test names(schema_run(-0.01u"eV"; obs=[:n_obs => cfg -> Float64(length(cfg))])) ==
              ["iter", "emax", "num_particles", "log_compression", "omega", "n_obs"]
    end

    @testset "_grand_potential closed form and the exact-zero addend" begin
        Random.seed!(86012)
        w = om_walker(3)
        w.energy = -0.5u"eV"
        # binary-exact values: -0.5 - (-0.25 * 3) = 0.25
        @test FreeBird.SamplingSchemes._grand_potential(w, -0.25u"eV") == 0.25u"eV"
        @test FreeBird.SamplingSchemes._grand_potential(w, 0.25u"eV") == -1.25u"eV"
        @test FreeBird.SamplingSchemes._grand_potential(w, 0.0u"eV") == w.energy
        # the sign of a negative zero survives the default's exact-zero addend
        w.energy = -0.0u"eV"
        @test signbit(ustrip(u"eV", FreeBird.SamplingSchemes._grand_potential(w, 0.0u"eV")))
        # the two-argument tie block agrees with the one-argument one at mu = 0 and
        # splits an energy plateau by particle number at mu != 0
        ws = [om_walker(n) for n in (2, 2, 3, 3)]
        for x in ws
            x.energy = 0.0u"eV"
        end
        @test FreeBird.SamplingSchemes._tie_block_length(ws) == 4
        @test FreeBird.SamplingSchemes._tie_block_length(ws, 0.0u"eV") == 4
        sort!(ws, by = x -> FreeBird.SamplingSchemes._grand_potential(x, -0.1u"eV"), rev=true)
        @test FreeBird.SamplingSchemes._tie_block_length(ws, -0.1u"eV") == 2
        @test ws[1].list_num_par[1] == 3
    end

    @testset "kernel: the Ω ceiling is a per-sector energy ceiling (eps = 0)" begin
        function three_particle_walker(mk, species, box)
            w = om_walker(0; mk=mk)
            for pos in ([2.0, 2.0, 2.0], [7.0, 7.0, 7.0], [2.0, 7.0, 2.0])
                FreeBird.AbstractWalkers.insert_particle!(w, SVector{3}(pos)u"Å", species)
            end
            w.energy = 0.0u"eV"
            return w
        end
        # Ensembles of short walks from N = 3 locate the boundary exactly: under
        # mu = -1 eV and a 3.5 eV ceiling the admissible sectors are N <= 3 (an
        # indicator evaluated on n instead of n + 1 would admit N = 4), under
        # mu = +1 eV and a -2.5 eV ceiling they are N >= 3 (an indicator on n
        # instead of n - 1 would admit N = 2); the mu = 0 controls show the same
        # ceilings admit both crossings when the sector shift is absent
        function ensemble(mu, ceiling; surface=false, n_walkers=200, steps=50)
            finals = Int[]
            ins = 0
            del = 0
            for _ in 1:n_walkers
                if surface
                    w = three_particle_walker(om_smkempty, :H, (10.0, 10.0, 15.0))
                    _, _, w, st = MC_grand_canonical_walk!(steps, w, om_cps0, ceiling, om_mksurf();
                        z0V=4.0, species=:H, p_move=0.4, p_insert=0.3, mu=mu)
                else
                    w = three_particle_walker(om_mkempty, :Ar, (12.0, 12.0, 12.0))
                    _, _, w, st = MC_grand_canonical_walk!(steps, w, om_lj0, ceiling;
                        z0V=4.0, species=:Ar, p_move=0.4, p_insert=0.3, step_size=0.8, mu=mu)
                end
                push!(finals, w.list_num_par[1])
                ins += st.insert_accepted
                del += st.delete_accepted
            end
            return finals, ins, del
        end
        Random.seed!(86020)
        f, ins, del = ensemble(-1.0u"eV", 3.5u"eV")
        @test maximum(f) == 3
        @test minimum(f) < 3
        @test ins > 0
        @test del > 0
        Random.seed!(86020)
        f0, _, _ = ensemble(0.0u"eV", 3.5u"eV")
        @test maximum(f0) > 3
        Random.seed!(86021)
        f, ins, del = ensemble(1.0u"eV", -2.5u"eV")
        @test minimum(f) == 3
        @test maximum(f) > 3
        @test ins > 0
        @test del > 0
        Random.seed!(86021)
        f0, _, _ = ensemble(0.0u"eV", 1.0e9u"eV")
        @test minimum(f0) < 3
        # the same four contracts through the surface-aware method
        Random.seed!(86022)
        f, _, _ = ensemble(-1.0u"eV", 3.5u"eV"; surface=true)
        @test maximum(f) == 3
        @test minimum(f) < 3
        Random.seed!(86022)
        f0, _, _ = ensemble(0.0u"eV", 3.5u"eV"; surface=true)
        @test maximum(f0) > 3
        Random.seed!(86023)
        f, _, _ = ensemble(1.0u"eV", -2.5u"eV"; surface=true)
        @test minimum(f) == 3
        @test maximum(f) > 3
        Random.seed!(86023)
        f0, _, _ = ensemble(0.0u"eV", 1.0e9u"eV"; surface=true)
        @test minimum(f0) < 3
        # an explicit mu = 0 is the shipped kernel draw for draw
        function kernel_pair(mu_kw)
            Random.seed!(86024)
            wk = three_particle_walker(om_mkempty, :Ar, (12.0, 12.0, 12.0))
            r = MC_grand_canonical_walk!(500, wk, om_lj0, 1.0u"eV";
                z0V=4.0, species=:Ar, p_move=0.4, p_insert=0.3, step_size=0.8, mu_kw...)
            return r[1], r[2], wk.list_num_par[1], wk.configuration.position, r[4]
        end
        a = kernel_pair(NamedTuple())
        b = kernel_pair((mu=0.0u"eV",))
        @test a[1] === b[1]
        @test a[2] == b[2]
        @test a[3] == b[3]
        @test a[4] == b[4]
        @test a[5] == b[5]
        # with a ceiling that never binds, a nonzero mu changes nothing either
        function loose_pair(mu)
            Random.seed!(86025)
            wk = three_particle_walker(om_mkempty, :Ar, (12.0, 12.0, 12.0))
            r = MC_grand_canonical_walk!(500, wk, om_lj0, 1.0e9u"eV";
                z0V=4.0, species=:Ar, p_move=0.4, p_insert=0.3, step_size=0.8, mu=mu)
            return wk.list_num_par[1], wk.configuration.position, r[4]
        end
        c = loose_pair(0.0u"eV")
        d = loose_pair(-1.0u"eV")
        @test c[1] == d[1]
        @test c[2] == d[2]
        @test c[3] == d[3]
    end

    @testset "kernel stationarity under the Ω ceiling (seeded, calibrated)" begin
        # ideal gas at z0V = 3, mu = -1 eV, ceiling 3.5 eV: the reference measure
        # conditioned on N <= 3, masses proportional to (1, 3, 4.5, 4.5)
        Random.seed!(86001)
        ns = Int[]
        for _ in 1:150
            w = om_walker(0)
            w.energy = 0.0u"eV"
            MC_grand_canonical_walk!(400, w, om_lj0, 3.5u"eV"; z0V=3.0, species=:Ar, mu=-1.0u"eV")
            push!(ns, w.list_num_par[1])
        end
        @test all(0 .<= ns .<= 3)
        @test abs(mean(ns) - 51 / 26) < 0.39
        @test abs(var(ns) - (123 / 26 - (51 / 26)^2)) < 0.42
        # mu = +1 eV, ceiling -2.5 eV, n_max = 6: conditioned on 3 <= N <= 6,
        # masses proportional to (4.5, 3.375, 2.025, 1.0125); walkers start at N = 3
        Random.seed!(86002)
        ns = Int[]
        for _ in 1:150
            w = om_walker(3)
            w.energy = 0.0u"eV"
            MC_grand_canonical_walk!(400, w, om_lj0, -2.5u"eV"; z0V=3.0, species=:Ar,
                                     n_max=6, mu=1.0u"eV")
            push!(ns, w.list_num_par[1])
        end
        q = [4.5, 3.375, 2.025, 1.0125] ./ 10.9125
        m1 = sum(q .* (3:6))
        v1 = sum(q .* (3:6) .^ 2) - m1^2
        @test all(3 .<= ns .<= 6)
        @test abs(mean(ns) - m1) < 0.53
        @test abs(var(ns) - v1) < 0.59
    end

    function om_run(seed, mu; K=16, z0V=6.0, mc_steps=60, n_steps=120, n_max=typemax(Int64),
                    tag="r", obs=nothing, callback=nothing, ls=nothing, initialize=true,
                    potential=om_lj, allowed_fail_count=100_000)
        Random.seed!(seed)
        lsr = ls === nothing ? GenericAtomWalkers([AtomWalker{1}(om_mkempty()) for _ in 1:K], potential) : ls
        p = AtomisticIGRefGCNSParameters(mc_steps=mc_steps, reference_activity=(z0V / om_V)u"Å^-3",
                                         species=:Ar, allowed_fail_count=allowed_fail_count,
                                         n_max=n_max, chemical_potential=mu)
        df, lso, pout = ideal_gas_referenced_nested_sampling(
            lsr, p, n_steps, MCAtomGrandCanonicalMoves(), om_save(tag);
            observables=obs, dead_point_callback=callback, initialize=initialize)
        om_clean(tag)
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]
        return df, live_e, live_n, lso, pout
    end

    @testset "driver stream identity: an explicit zero is the default" begin
        df_d, e_d, n_d, _, _ = om_run(52552, 0.0u"eV")
        Random.seed!(52552)
        ls = GenericAtomWalkers([AtomWalker{1}(om_mkempty()) for _ in 1:16], om_lj)
        p = AtomisticIGRefGCNSParameters(mc_steps=60, reference_activity=(6.0 / om_V)u"Å^-3",
                                         species=:Ar, allowed_fail_count=100_000)
        df_0, lso_0, _ = ideal_gas_referenced_nested_sampling(
            ls, p, 120, MCAtomGrandCanonicalMoves(), om_save("z"))
        om_clean("z")
        @test df_d == df_0
        @test e_d == [ustrip(u"eV", w.energy) for w in lso_0.walkers]
        @test n_d == [w.list_num_par[1] for w in lso_0.walkers]
        # the surface path likewise
        function sf_pair(mu_kw)
            Random.seed!(84004)
            lss = om_sliveset(8)
            ps = AtomisticIGRefGCNSParameters(; mc_steps=40, reference_activity=(4.0 / om_sV)u"Å^-3",
                                              species=:H, allowed_fail_count=1000, n_max=6, mu_kw...)
            dfs, lsos, _ = ideal_gas_referenced_nested_sampling(
                lss, ps, 60, MCAtomGrandCanonicalMoves(), om_save("s"))
            om_clean("s")
            return dfs, [ustrip(u"eV", w.energy) for w in lsos.walkers]
        end
        sa = sf_pair(NamedTuple())
        sb = sf_pair((chemical_potential=0.0u"eV",))
        @test sa[1] == sb[1]
        @test sa[2] == sb[2]
    end

    @testset "ordering contract at mu != 0 on the dilute interacting gas (seed 86030)" begin
        # mu = -0.002 eV: bound pairs (E about -0.01 eV) sit below the exact
        # N = 1 plateau at Ω = 0.002 eV, so the descent has an interior to refill
        # from and every step records a row
        mu = -0.002u"eV"
        df1, e1, n1, _, _ = om_run(86030, mu)
        df2, e2, n2, _, _ = om_run(86030, mu)
        @test df1 == df2
        @test e1 == e2
        @test n1 == n2
        @test nrow(df1) >= 60
        @test length(unique(df1.iter)) == nrow(df1)
        @test issorted(df1.omega, rev=true)
        # the energy column is not the ordering scalar any more
        @test !issorted(df1.emax, rev=true)
        # the recorded Ω is the step's own expression on the recorded (E, N)
        @test df1.omega == df1.emax .- ustrip(u"eV", mu) .* df1.num_particles
        @test all(df1.log_compression .< 0.0)
    end

    @testset "observable pre-sort and dead-point pairing at mu != 0 (seed 86031)" begin
        seen = Int[]
        df, _, _, _, _ = om_run(86031, -0.002u"eV"; K=8, z0V=4.0, mc_steps=40, n_steps=40,
                                obs=[:n_obs => cfg -> Float64(length(cfg))],
                                callback=(iter, walker) -> push!(seen, walker.list_num_par[1]))
        @test names(df) == ["iter", "emax", "num_particles", "log_compression", "omega", "n_obs"]
        @test df.n_obs == Float64.(df.num_particles)
        @test seen == df.num_particles
        @test issorted(df.omega, rev=true)
    end

    @testset "atoms of the ordering scalar end the dilute descent in a stall, never a throw (seed 86032)" begin
        # mu = -0.02 eV on the dilute gas at K = 8: the exact N = 1 level
        # (Ω = 0.02 eV, every isolated particle) and the empty configuration
        # (Ω = 0) are atoms of Ω. When the N = 1 level is the top block, its
        # survivors are empties, whose clones stay below it only by accepting
        # no move, so the refill exhausts its budget and the live set shrinks;
        # when either atom holds the whole live set, the step attempts no walk.
        # Either way the driver stalls on an all-tied live set and never draws
        # a parent from an empty range (the shipped code threw there).
        df, live_e, live_n, lso, pout = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            om_run(86032, -0.02u"eV"; K=8, z0V=4.0, mc_steps=40, n_steps=400, allowed_fail_count=20)
        end
        @test pout.fail_count == 20
        @test 1 <= length(lso.walkers) <= 8
        omegas = [FreeBird.SamplingSchemes._grand_potential(w, -0.02u"eV") for w in lso.walkers]
        @test all(==(omegas[1]), omegas)
        @test all(n -> n <= 1, live_n)
        @test issorted(df.omega, rev=true)
        @test df.omega == df.emax .- (-0.02) .* df.num_particles
    end

    @testset "ideal-gas descent under Ω ordering closes against the Poisson closed forms (seed 86003)" begin
        # U = 0 everywhere, so Ω = -μN: at mu < 0 the descent evicts the particle-number
        # plateaus from the top down and ends on K empties (the terminal atom); the
        # compression ledger then encodes the reference law's N-marginal
        z0V = 4.0
        z0 = (z0V / om_V)u"Å^-3"
        T = 300.0
        mu = -0.02u"eV"
        df, live_e, live_n, _, _ = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            om_run(86003, mu; K=256, z0V=z0V, mc_steps=20, n_steps=4000, potential=IdealGasParameters(),
                   allowed_fail_count=10, tag="ig")
        end
        @test nrow(df) > 256
        @test all(df.emax .== 0.0)
        @test issorted(df.omega, rev=true)
        @test issorted(df.num_particles, rev=true)
        @test df.omega == -ustrip(u"eV", mu) .* df.num_particles
        @test all(iszero, live_n)
        @test all(iszero, live_e)
        mus = [om_mu_for(0.7 * z0V, T), om_mu_for(z0V, T), om_mu_for(1.3 * z0V, T)]
        st = gc_thermodynamic_stats_ideal_ref(df, om_Vq, om_mass, z0, mus, [T]u"K";
                                              live_emax=live_e, live_numbers=live_n)
        # at the reference activity the shells and the tail telescope to the full
        # reference mass, whatever ordering produced them
        @test abs(st.logXi[2, 1] - z0V) < 1e-10
        @test all(st.mean_U .== 0.0)
        # the Ω-ordered compression reproduces the Poisson N-marginal: activity
        # linearity and the occupancy moments at every grid point (calibrated)
        @test abs(st.logXi[1, 1] - 0.7 * z0V) < 0.24
        @test abs(st.logXi[3, 1] - 1.3 * z0V) < 0.24
        for (i, zV) in enumerate((0.7 * z0V, z0V, 1.3 * z0V))
            @test abs(st.mean_N[i, 1] - zV) < 1.25
            @test abs(st.var_N[i, 1] - zV) < 4.8
            pois = [exp(n * log(zV) - zV - om_logfact(n)) for n in st.N_support]
            @test 0.5 * sum(abs.(st.p_N[i, 1, :] .- pois)) < 0.28
        end
        # the bounded construction: the same descent on the conditional Poisson
        n_cap = 6
        dfb, live_eb, live_nb, _, _ = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            om_run(86004, mu; K=256, z0V=z0V, mc_steps=20, n_steps=4000, potential=IdealGasParameters(),
                   allowed_fail_count=10, n_max=n_cap, tag="igb")
        end
        @test maximum(dfb.num_particles) <= n_cap
        @test issorted(dfb.num_particles, rev=true)
        @test all(iszero, live_nb)
        stb = gc_thermodynamic_stats_ideal_ref(dfb, om_Vq, om_mass, z0, mus, [T]u"K";
                                               live_emax=live_eb, live_numbers=live_nb, n_max=n_cap)
        @test abs(stb.logXi[2, 1] - om_lpps(z0V, n_cap)) < 1e-10
        @test abs(stb.logXi[1, 1] - om_lpps(0.7 * z0V, n_cap)) < 0.24
        @test abs(stb.logXi[3, 1] - om_lpps(1.3 * z0V, n_cap)) < 0.24
        # mu > 0 under a cap: the descent climbs the plateaus and ends on the cap
        n_cap_p = 8
        dfp, live_ep, live_np, _, _ = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            om_run(86005, -mu; K=256, z0V=z0V, mc_steps=20, n_steps=4000, potential=IdealGasParameters(),
                   allowed_fail_count=10, n_max=n_cap_p, tag="igp")
        end
        @test all(dfp.emax .== 0.0)
        @test issorted(dfp.omega, rev=true)
        @test issorted(dfp.num_particles)
        @test all(==(n_cap_p), live_np)
        stp = gc_thermodynamic_stats_ideal_ref(dfp, om_Vq, om_mass, z0, mus, [T]u"K";
                                               live_emax=live_ep, live_numbers=live_np, n_max=n_cap_p)
        @test abs(stp.logXi[2, 1] - om_lpps(z0V, n_cap_p)) < 1e-10
        @test abs(stp.logXi[1, 1] - om_lpps(0.7 * z0V, n_cap_p)) < 0.24
        @test abs(stp.logXi[3, 1] - om_lpps(1.3 * z0V, n_cap_p)) < 0.24
    end

    @testset "terminal-atom stall over a frozen substrate (seed 86040)" begin
        # |mu| far above every interaction energy: the Ω order is the particle-number
        # order, the descent empties every walker, and the empty configuration is an
        # atom of Ω (every empty walker carries exactly energy_frozen_part)
        Random.seed!(86040)
        walkers = AtomWalker{1}[]
        for n in (0, 1, 1, 2, 2, 3, 3, 4)
            push!(walkers, om_walker(n; species=:H, box=(10.0, 10.0, 15.0), mk=om_smkempty))
        end
        ls = LJSurfaceWalkers(walkers, om_cps, om_mksurf(); assign_energy=true)
        e_frozen = ls.surface.energy_frozen_part
        mu = -1.0u"eV"
        p = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / om_sV)u"Å^-3",
                                         species=:H, allowed_fail_count=5, chemical_potential=mu)
        df, lso, pout = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ideal_gas_referenced_nested_sampling(ls, p, 500, MCAtomGrandCanonicalMoves(), om_save("atom");
                                                 initialize=false)
        end
        om_clean("atom")
        @test pout.fail_count == 5
        @test nrow(df) >= 7
        @test all(w.list_num_par[1] == 0 for w in lso.walkers)
        @test all(w.energy == e_frozen for w in lso.walkers)
        @test issorted(df.omega, rev=true)
        @test issorted(df.num_particles, rev=true)
        @test df.omega == df.emax .- ustrip(u"eV", mu) .* df.num_particles
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]
        st = gc_thermodynamic_stats_ideal_ref(df, 1500.0u"Å^3", 1.008u"u", (4.0 / om_sV)u"Å^-3",
                                              [mu], [300.0]u"K"; live_emax=live_e, live_numbers=live_n)
        @test st.mean_N[1, 1] < 1e-10
        @test abs(st.mean_U[1, 1] - ustrip(u"eV", e_frozen)) < 1e-12
        # the vacuous effective sample size of an atom-dominated target: the K
        # identical tail walkers share the weight
        @test abs(st.N_eff[1, 1] - 8.0) < 1e-6
        # a continuation of the stalled set stalls again without a recorded row
        df2, lso2, _ = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ideal_gas_referenced_nested_sampling(lso, p, 50, MCAtomGrandCanonicalMoves(), om_save("atom2");
                                                 initialize=false)
        end
        om_clean("atom2")
        @test nrow(df2) == 0
        @test all(w.list_num_par[1] == 0 for w in lso2.walkers)
    end

    @testset "an all-tied live set is terminal under a nonzero mu (seed 86060)" begin
        # K identical empties over the substrate sit on one Ω; bound single
        # adsorbates lie below it at this mu (E_1 is about -0.01 eV against
        # mu = -0.001 eV) and are easy to find from an empty parent. At mu = 0 the
        # shipped path culls through the atom and records rows; under the
        # nonzero mu the step attempts no walk and the stall fires with the
        # live set intact, its mass charged to the atom.
        function tied_run(mu)
            Random.seed!(86060)
            ls = om_sliveset(6)
            p = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / om_sV)u"Å^-3",
                                             species=:H, allowed_fail_count=8,
                                             chemical_potential=mu)
            df, lso, pout = ideal_gas_referenced_nested_sampling(
                ls, p, 30, MCAtomGrandCanonicalMoves(), om_save("tied"); initialize=false)
            om_clean("tied")
            return df, lso, pout
        end
        df0, lso0, _ = tied_run(0.0u"eV")
        @test nrow(df0) > 0
        @test any(w.list_num_par[1] > 0 for w in lso0.walkers)
        dfm, lsom, poutm = @test_logs (:warn, r"stop_on_stall") match_mode=:any tied_run(-0.001u"eV")
        @test nrow(dfm) == 0
        @test poutm.fail_count == 8
        @test all(w.list_num_par[1] == 0 for w in lsom.walkers)
        @test length(lsom.walkers) == 6
        # the same live set with one walker below the atom is an ordinary
        # plateau: the evictions and the refill proceed and rows are recorded
        Random.seed!(86061)
        ls1 = om_sliveset(6)
        # one adsorbate at the pair minimum 2^(1/6) sigma = 2.81 A above a
        # surface atom (the other three sit 5 A away, beyond the cutoff)
        w1 = AtomWalker{1}(om_smkempty())
        FreeBird.AbstractWalkers.insert_particle!(w1, SVector(2.5, 2.5, 2.0 + 2.0^(1 / 6) * 2.5)u"Å", :H)
        ls1.walkers[1] = w1
        FreeBird.SamplingSchemes._reanchor_igref_energies!(ls1)
        p1 = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / om_sV)u"Å^-3",
                                          species=:H, allowed_fail_count=8,
                                          chemical_potential=-0.001u"eV")
        @test FreeBird.SamplingSchemes._grand_potential(ls1.walkers[1], -0.001u"eV") < 0.0u"eV"
        df1, lso1, _ = ideal_gas_referenced_nested_sampling(
            ls1, p1, 30, MCAtomGrandCanonicalMoves(), om_save("tied1"); initialize=false)
        om_clean("tied1")
        @test nrow(df1) >= 5
        @test count(==(0), df1.num_particles) == 5
    end

    @testset "a continuation carries the plateau refill target across a chunk boundary (seed 86070)" begin
        # K = 8 over the substrate: four empties on the atom Ω = 0 (the top block
        # under mu = -0.001 eV) above four bound adsorbates with distinct heights.
        # A chunked driver that re-enters mid-eviction with the parameters it
        # carried keeps the target 8 and refills to 8 once the plateau is
        # exhausted; a re-entry with fresh parameters re-derives the target from
        # the reduced live count and never refills the walkers evicted before
        # the boundary
        function mid_eviction_liveset()
            ls = om_sliveset(8)
            for (i, z) in enumerate((4.80, 4.86, 4.92, 4.98))
                w = AtomWalker{1}(om_smkempty())
                FreeBird.AbstractWalkers.insert_particle!(w, SVector(2.5, 2.5, z)u"Å", :H)
                ls.walkers[4 + i] = w
            end
            FreeBird.SamplingSchemes._reanchor_igref_energies!(ls)
            return ls
        end
        mu = -0.001u"eV"
        mkp() = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / om_sV)u"Å^-3",
                                             species=:H, allowed_fail_count=200,
                                             chemical_potential=mu)
        Random.seed!(86070)
        ls = mid_eviction_liveset()
        p = mkp()
        df1, lso1, p1 = ideal_gas_referenced_nested_sampling(ls, p, 2, MCAtomGrandCanonicalMoves(), om_save("pt1");
                                                             initialize=false)
        om_clean("pt1")
        @test nrow(df1) == 2
        @test all(==(0), df1.num_particles)
        @test length(lso1.walkers) == 6
        @test p1.plateau_refill_target == 8
        # carried parameters: the two remaining empties are evicted, then the
        # refill restores the live set to 8
        df2, lso2, p2 = ideal_gas_referenced_nested_sampling(lso1, p1, 6, MCAtomGrandCanonicalMoves(), om_save("pt2");
                                                             initialize=false)
        om_clean("pt2")
        @test length(lso2.walkers) == 8
        @test p2.plateau_refill_target == 0
        # fresh parameters at the same boundary: the target is re-derived as 6
        Random.seed!(86070)
        lsb = mid_eviction_liveset()
        dfb1, lsob1, _ = ideal_gas_referenced_nested_sampling(lsb, mkp(), 2, MCAtomGrandCanonicalMoves(), om_save("pt3");
                                                              initialize=false)
        om_clean("pt3")
        dfb2, lsob2, _ = ideal_gas_referenced_nested_sampling(lsob1, mkp(), 6, MCAtomGrandCanonicalMoves(), om_save("pt4");
                                                              initialize=false)
        om_clean("pt4")
        @test length(lsob2.walkers) == 6
    end

    @testset "block continuation carries mu (seed 86050)" begin
        Random.seed!(86050)
        ls = om_sliveset(8)
        mu = -0.005u"eV"
        p = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(4.0 / om_sV)u"Å^-3",
                                         species=:H, allowed_fail_count=1000, n_max=6,
                                         chemical_potential=mu)
        df1, lso1, _ = ideal_gas_referenced_nested_sampling(ls, p, 60, MCAtomGrandCanonicalMoves(), om_save("c1"))
        om_clean("c1")
        df2, lso2, _ = ideal_gas_referenced_nested_sampling(lso1, p, 40, MCAtomGrandCanonicalMoves(), om_save("c2");
                                                            initialize=false)
        om_clean("c2")
        @test nrow(df2) > 0
        @test df2.iter[1] == df1.iter[end] + 1
        @test df2.omega[1] <= df1.omega[end]
        @test issorted(vcat(df1.omega, df2.omega), rev=true)
        @test df2.omega == df2.emax .- ustrip(u"eV", mu) .* df2.num_particles
    end

    @testset "reference_activity_temperature" begin
        for (z0V, Vc) in ((16.0, 2091.124), (4.0, 2091.124), (3.0, 1000.0))
            z0 = (z0V / Vc)u"Å^-3"
            TL = reference_activity_temperature(z0, om_mass)
            @test unit(TL) == u"K"
            lam3 = FreeBird.AnalysisTools._thermal_wavelength(om_mass, TL)^3
            @test abs(ustrip(Unitful.NoUnits, z0 * lam3) - 1.0) < 1e-12
            # Λ^3 scales as T^(-3/2): eight times the activity quadruples the temperature
            @test isapprox(reference_activity_temperature(8 * z0, om_mass), 4 * TL; rtol=1e-12)
        end
        @test_throws ArgumentError reference_activity_temperature(0.0u"Å^-3", om_mass)
        @test_throws ArgumentError reference_activity_temperature(-1.0u"Å^-3", om_mass)
    end
end
