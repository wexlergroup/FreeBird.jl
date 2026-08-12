# RNG-neutrality and closed-form regression coverage for the continuous-space
# grand-canonical route.
#
# Calibration ledger:
# - (a) The canonical-ladder pins were recorded by running the identical
#   fixture on dev @ 6b9d030 and on that change set's bytes: the two streams
#   agree digit-for-digit on one machine (rows 60, last emax
#   -0.024149711142342864 eV, live-energy sum -0.20617927437954883 eV, last
#   log-compression -0.11778303565638351). Across architectures and compiler
#   versions the same trajectory reproduces these values only to the last one
#   or two ulps (ubuntu-x64 under Julia 1.12 measured -0.024149711142342878
#   against the macOS-aarch64 recording), so the pins gate at rtol 1e-12: a
#   genuine stream change produces order-one deviations and still fails
#   loudly, while compiler-level rounding drift along the same trajectory
#   passes.
# - (b) Ideal-gas pipeline (K = 512, z0V = 4, 50 requested steps stalling at 5,
#   mu grid at target zV = {2.8, 4, 5.2}, seeds {1,2,3,81811}): per-point max
#   devs logXi 0.063, mean_N {0.084, 0.139, 0.386}, var_N {0.142, 0.498, 1.518};
#   gates ship at >= 3x per point (logXi 0.2; mean_N 0.25, 0.45, 1.2; var_N
#   0.45, 1.5, 4.6). Observed N_eff floor 361 behind the > 100 assertion. At
#   the reference activity the assembled logXi equals z0V exactly (the
#   reference-mass prefactor and the tail normalization cancel bitwise);
#   observed deviation 0.0 on every seed.
# - (c) Kernel energy drift over 2e4 insert/delete/move steps on a dilute
#   interacting gas (seeds {1,2,3,81812}): max 3.1e-11 eV; bound ships at 5e-10.
# - (d) Cap-inactivity run (K = 16, z0V = 6, 120 steps, dilute interacting gas,
#   seeds {1,2,3,81813}): observed max particle count 15 (ledger and live set
#   alike) against the 50 bound, a 3x-plus margin on the support edge.
# - (e) The save layer cannot currently persist a zero-atom frame (the ExtXYZ
#   writer fails converting an empty species table); the testset pins that
#   failure as the documented disclosure, so a future ExtXYZ that learns to
#   write empty frames surfaces here and the persistence testset can be added.
@testset "atomistic grand-canonical regression coverage" begin
    using Random

    box12 = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
    box10 = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    pbc = (true, true, true)
    seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box12, pbc))
    mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                           empty(position(seed_at, :)), empty(species(seed_at, :)),
                           empty(mass(seed_at, :)))
    V = 1728.0
    Vq = 1728.0u"Å^3"
    mass_ar = 39.948u"u"
    kb = 8.617333262e-5
    lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass_ar, T * u"K"))
    mu_for(zV, T) = (kb * T * (log(zV / V) + 3 * log(lam(T)))) * u"eV"
    mksave(tag; n_snap=10^7) = SaveEveryN(df_filename="_gcreg_$(tag).csv",
                                          wk_filename="_gcreg_$(tag).traj.extxyz",
                                          ls_filename="_gcreg_$(tag).ls.extxyz",
                                          n_traj=10^7, n_snap=n_snap, n_info=10^7)
    clean(tag) = for f in ["_gcreg_$(tag).csv", "_gcreg_$(tag).traj.extxyz", "_gcreg_$(tag).ls.extxyz"]
        rm(f, force=true)
    end

    @testset "RNG-stream neutrality: canonical-ladder digit pins" begin
        Random.seed!(71711)
        walkers = AtomWalker{1}[]
        for _ in 1:8
            coords = [[rand() * 10.0, rand() * 10.0, rand() * 10.0] for _ in 1:3]
            push!(walkers, AtomWalker{1}(FastSystem(atomic_system(
                [:Ar => SVector{3}(c)u"Å" for c in coords], box10, pbc))))
        end
        ls = LJAtomWalkers(walkers, LJParameters(epsilon=0.01, sigma=2.5))
        p = NestedSamplingParameters(mc_steps=40, initial_step_size=0.5, step_size=0.5,
                                     step_size_lo=0.01, step_size_up=2.0,
                                     allowed_fail_count=1000)
        save = mksave("ab")
        df, lso, _ = nested_sampling(ls, p, 60, MCRandomWalkClone(), save)
        clean("ab")
        @test nrow(df) == 60
        # rtol 1e-12: tolerate compiler-level rounding drift along the same
        # trajectory (see the calibration ledger); a stream change fails loudly
        @test isapprox(df.emax[end], -0.024149711142342864; rtol=1e-12)
        @test isapprox(sum(ustrip(u"eV", w.energy) for w in lso.walkers), -0.20617927437954883; rtol=1e-12)
        @test isapprox(df.log_compression[end], -0.11778303565638351; rtol=1e-12)
    end

    @testset "end-to-end ideal-gas pipeline closes against the closed forms (seed 81811)" begin
        z0V = 4.0
        z0 = (z0V / V)u"Å^-3"
        T = 300.0
        Random.seed!(81811)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:512], IdealGasParameters())
        params = AtomisticIGRefGCNSParameters(mc_steps=10, reference_activity=z0,
                                              species=:Ar, allowed_fail_count=5)
        save = mksave("pipe")
        df, lso, _ = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ideal_gas_referenced_nested_sampling(ls, params, 50, MCAtomGrandCanonicalMoves(), save)
        end
        clean("pipe")
        @test nrow(df) == 0
        @test names(df) == ["iter", "emax", "num_particles", "log_compression"]
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]
        @test all(iszero, live_e)
        mus = [mu_for(0.7 * z0V, T), mu_for(z0V, T), mu_for(1.3 * z0V, T)]
        st = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                              live_emax=live_e, live_numbers=live_n)
        # at the reference activity the reference-mass prefactor and the tail
        # normalization cancel, so the assembled logXi equals z0V up to rounding
        @test abs(st.logXi[2, 1] - z0V) < 1e-12
        # the activity-linearity identity logXi(z) = zV, one gate testing the
        # prefactor, the reweighting factor, and the 1/N! bookkeeping together
        @test abs(st.logXi[1, 1] - 0.7 * z0V) < 0.2
        @test abs(st.logXi[3, 1] - 1.3 * z0V) < 0.2
        # occupation mean and variance close on every grid point (per-point gates)
        @test abs(st.mean_N[1, 1] - 0.7 * z0V) < 0.25
        @test abs(st.mean_N[2, 1] - z0V) < 0.45
        @test abs(st.mean_N[3, 1] - 1.3 * z0V) < 1.2
        @test abs(st.var_N[1, 1] - 0.7 * z0V) < 0.45
        @test abs(st.var_N[2, 1] - z0V) < 1.5
        @test abs(st.var_N[3, 1] - 1.3 * z0V) < 4.6
        @test all(st.mean_U .== 0.0)
        @test minimum(st.N_eff) > 100.0
    end

    @testset "kernel incremental-energy drift pin (seed 81812)" begin
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
        Random.seed!(81812)
        w = AtomWalker{1}(mkempty())
        w.energy = 0.0u"eV"
        MC_grand_canonical_walk!(20000, w, lj, 1.0e5u"eV"; z0V=4.0, species=:Ar, step_size=1.0)
        recomputed = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen)
        @test abs(ustrip(u"eV", w.energy - recomputed)) < 5.0e-10
        @test length(w.configuration) == w.list_num_par[1]
    end

    @testset "particle-number cap stays structurally inactive" begin
        # the sampler's parameter struct carries no cap at all: a cap would
        # silently truncate the reference measure's unbounded support
        @test !hasfield(AtomisticIGRefGCNSParameters, :n_max)
        Random.seed!(81813)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:16],
                                LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5))
        params = AtomisticIGRefGCNSParameters(mc_steps=60, reference_activity=(6.0 / V)u"Å^-3",
                                              species=:Ar, allowed_fail_count=50)
        save = mksave("cap")
        df, lso, _ = ideal_gas_referenced_nested_sampling(ls, params, 120,
                                                          MCAtomGrandCanonicalMoves(), save)
        clean("cap")
        @test maximum(df.num_particles) < 50
        @test maximum(w.list_num_par[1] for w in lso.walkers) < 50
    end

    @testset "empty-walker persistence: documented save-layer limitation" begin
        ws = [AtomWalker{1}(mkempty()),
              AtomWalker{1}(FastSystem(atomic_system([:Ar => [2.0, 2.0, 2.0]u"Å"], box12, pbc)))]
        ls = GenericAtomWalkers(ws, LJParameters(epsilon=0.0))
        save = mksave("io"; n_snap=1)
        # DISCLOSURE: the ExtXYZ-backed save layer cannot represent a zero-atom
        # frame today (species-table conversion fails). This pin keeps the
        # limitation documented; when it stops throwing, replace it with a
        # write-and-read-back persistence testset.
        @test_throws MethodError FreeBird.FreeBirdIO.write_ls_every_n(ls, 1, save)
        clean("io")
    end
end
