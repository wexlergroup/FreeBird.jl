# Surface-aware ideal-gas-referenced grand-canonical route (LJSurfaceWalkers).
#
# Calibration ledger:
# - Athermal stationary law (kernel under an effectively infinite ceiling,
#   150 independent walkers x 400 steps, z0V = 3, seeds {1,2,3,84003}): the
#   stationary law is the reference measure itself, independent of the
#   potential; max devs mean 0.140, var 0.343 (both at the shipped seed);
#   gates ship at >= 3x (0.42 and 1.05).
# - Capstone bounded cross-route (surface igref, z0V = 4, n_max = 4, K = 16,
#   200 steps vs fixed-N stitching of canonical LJSurfaceWalkers ladders over
#   N = 0..4, K = 16, 120 steps; per-T mu grids solved from target zV in
#   {2, 4}; igref seeds {1,2,3,84621} against fixed-N seeds
#   {1001,1002,1003,85621}): max devs logXi 0.333, mean_N 0.720 (both on the seed-3 pair),
#   N_eff floor observed 11.8 (binding on the seed-1 calibration pair);
#   gates ship at >= 3x (1.0 and 2.2) with an N_eff floor at 3.5. During
#   calibration the
#   gate centers were cross-checked against a brute-force uniform-sampling
#   oracle for the truncated Xi (not shipped).
# - The transparent-surface weld is argued exact in-place: with eps_as = 0
#   every surface pair term is exactly 0.0 eV (4*0*(finite) with no overflow at
#   physical distances), the adsorbate-adsorbate partial sums iterate in
#   identical order in both kernel methods, and the two methods consume
#   identical RNG streams, so the trajectories are bitwise equal; the rtol
#   1e-12 fallback is disclosed should an x64 leg falsify the argument.
@testset "surface-aware ideal-gas-referenced route" begin
    using Random

    sf_box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
    sf_pbc = (true, true, false)
    sf_V = 1500.0
    sf_surf_sys = FastSystem(atomic_system(
        [:H => [2.5, 2.5, 2.0]u"Å", :H => [7.5, 2.5, 2.0]u"Å",
         :H => [2.5, 7.5, 2.0]u"Å", :H => [7.5, 7.5, 2.0]u"Å"], sf_box, sf_pbc))
    sf_mksurf() = AtomWalker(deepcopy(sf_surf_sys); freeze_species=[:H])
    sf_mkempty() = FastSystem(cell_vectors(sf_surf_sys), periodicity(sf_surf_sys),
                              empty(position(sf_surf_sys, :)), empty(species(sf_surf_sys, :)),
                              empty(mass(sf_surf_sys, :)))
    # cutoff 1.8 sigma = 4.5 A stays inside every half cell edge: no
    # minimum-image warning noise in this file
    sf_lj_aa = LJParameters(epsilon=0.001, sigma=2.5, cutoff=1.8, shift=true)
    sf_lj_as = LJParameters(epsilon=0.003, sigma=2.5, cutoff=1.8, shift=true)
    sf_lj_ss = LJParameters(epsilon=0.01, sigma=2.5, cutoff=1.8, shift=true)
    sf_cps = CompositeParameterSets(2, [sf_lj_aa, sf_lj_as, sf_lj_ss])
    sf_save(tag) = SaveEveryN(df_filename="_igsurf_$(tag).csv",
                              wk_filename="_igsurf_$(tag).traj.extxyz",
                              ls_filename="_igsurf_$(tag).ls.extxyz",
                              n_traj=10^7, n_snap=10^7, n_info=10^7)
    sf_clean(tag) = for f in ["_igsurf_$(tag).csv", "_igsurf_$(tag).traj.extxyz",
                              "_igsurf_$(tag).ls.extxyz"]
        rm(f, force=true)
    end
    function sf_liveset(K)
        LJSurfaceWalkers([AtomWalker{1}(sf_mkempty()) for _ in 1:K], sf_cps, sf_mksurf();
                         assign_energy=true)
    end

    @testset "kernel validation and channel absence" begin
        surf = sf_mksurf()
        w = AtomWalker{1}(sf_mkempty())
        w.energy = 0.0u"eV"
        # a parameter set that is not two-component is rejected
        cps3 = CompositeParameterSets(3, [sf_lj_aa, sf_lj_aa, sf_lj_aa,
                                          sf_lj_aa, sf_lj_aa, sf_lj_aa])
        @test_throws ArgumentError MC_grand_canonical_walk!(5, w, cps3, 1.0u"eV", surf;
                                                            z0V=3.0, species=:H)
        # a frozen adsorbate walker is rejected
        wf = AtomWalker{1}(sf_mkempty())
        wf.frozen = [true]
        @test_throws ArgumentError MC_grand_canonical_walk!(5, wf, sf_cps, 1.0u"eV", surf;
                                                            z0V=3.0, species=:H)
        # a surface in a different cell is rejected
        other_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
        other_surf = AtomWalker(FastSystem(atomic_system(
            [:H => [2.5, 2.5, 2.0]u"Å"], other_box, sf_pbc)); freeze_species=[:H])
        @test_throws ArgumentError MC_grand_canonical_walk!(5, w, sf_cps, 1.0u"eV", other_surf;
                                                            z0V=3.0, species=:H)
        # the cavity-biased channel does not exist on the surface method: a
        # kwarg MethodError, never a silently ignored option
        @test_throws MethodError MC_grand_canonical_walk!(5, w, sf_cps, 1.0u"eV", surf;
                                                          z0V=3.0, species=:H, p_bias=0.4)
        @test_throws MethodError MC_grand_canonical_walk!(5, w, sf_cps, 1.0u"eV", surf;
                                                          z0V=3.0, species=:H,
                                                          bias_radius=2.0, bias_grid=8)
    end

    @testset "transparent-surface weld (eps_as = 0, seed 84001)" begin
        lj0as = LJParameters(epsilon=0.0, sigma=2.5, cutoff=1.8, shift=false)
        cps0 = CompositeParameterSets(2, [sf_lj_aa, lj0as, sf_lj_ss])
        surf = sf_mksurf()
        w1 = AtomWalker{1}(sf_mkempty())
        w1.energy = 0.0u"eV"
        w2 = AtomWalker{1}(sf_mkempty())
        w2.energy = 0.0u"eV"
        Random.seed!(84001)
        r1 = MC_grand_canonical_walk!(300, w1, cps0, 1.0u"eV", surf; z0V=3.0, species=:H)
        Random.seed!(84001)
        r2 = MC_grand_canonical_walk!(300, w2, sf_lj_aa, 1.0u"eV"; z0V=3.0, species=:H)
        @test w1.energy == w2.energy
        @test w1.list_num_par[1] == w2.list_num_par[1]
        @test w1.configuration.position == w2.configuration.position
        @test r1[4] == r2[4]
    end

    @testset "incremental-energy drift pin (seed 84002)" begin
        Random.seed!(84002)
        surf = sf_mksurf()
        w = AtomWalker{1}(sf_mkempty())
        w.energy = 0.0u"eV"
        MC_grand_canonical_walk!(20000, w, sf_cps, 1.0e5u"eV", surf;
                                 z0V=4.0, species=:H, step_size=1.0)
        recomputed = interacting_energy(w.configuration, sf_cps, w.list_num_par,
                                        w.frozen, surf.configuration)
        @test abs(ustrip(u"eV", w.energy - recomputed)) < 5.0e-10
        @test length(w.configuration) == w.list_num_par[1]
    end

    @testset "athermal stationary reference law (seeded, calibrated)" begin
        # Under an effectively infinite ceiling the kernel's stationary law is
        # the reference measure itself, whatever the potential: Poisson counts
        Random.seed!(84003)
        surf = sf_mksurf()
        ns = Int[]
        for _ in 1:150
            w = AtomWalker{1}(sf_mkempty())
            w.energy = 0.0u"eV"
            MC_grand_canonical_walk!(400, w, sf_cps, 1.0e9u"eV", surf;
                                     z0V=3.0, species=:H)
            push!(ns, w.list_num_par[1])
        end
        @test abs(mean(ns) - 3.0) < 0.42
        @test abs(var(ns) - 3.0) < 1.05
    end

    @testset "surface driver internals" begin
        ls = sf_liveset(4)
        params = AtomisticIGRefGCNSParameters(reference_activity=(3.0 / sf_V)u"Å^-3",
                                              species=:H)
        @test FreeBird.SamplingSchemes._atomistic_igref_z0V(ls, params) ≈ 3.0 atol = 1e-12
        # a frozen ADSORBATE component is rejected even on the surface path
        ls.walkers[1].frozen = [true]
        @test_throws ArgumentError FreeBird.SamplingSchemes._atomistic_igref_z0V(ls, params)
        ls.walkers[1].frozen = [false]
        # initialization inherits the surface's frozen-part energy and assigns
        # against the surface
        Random.seed!(84010)
        ls.surface.energy_frozen_part = 0.123u"eV"
        FreeBird.SamplingSchemes._init_atomistic_igref_walkers!(ls, params, 3.0)
        @test all(w.energy_frozen_part == 0.123u"eV" for w in ls.walkers)
        @test all(w.iter == 0 for w in ls.walkers)
        for w in ls.walkers
            fresh = interacting_energy(w.configuration, ls.potential, w.list_num_par,
                                       w.frozen, ls.surface.configuration) + w.energy_frozen_part
            @test w.energy == fresh
        end
        # driver-level surface checks: a not-fully-frozen surface and a
        # surface in a different cell are rejected by the validator
        lsn = sf_liveset(2)
        lsn.surface.frozen = [false]
        @test_throws ArgumentError FreeBird.SamplingSchemes._atomistic_igref_z0V(lsn, params)
        other_box2 = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
        far_surf = AtomWalker(FastSystem(atomic_system(
            [:H => [2.5, 2.5, 2.0]u"Å"], other_box2, sf_pbc)); freeze_species=[:H])
        lsm = LJSurfaceWalkers([AtomWalker{1}(sf_mkempty()) for _ in 1:2], sf_cps, far_surf;
                               assign_energy=false)
        @test_throws ArgumentError FreeBird.SamplingSchemes._atomistic_igref_z0V(lsm, params)
        # non-surface-aware channels are rejected at step entry
        ls2 = sf_liveset(4)
        params2 = AtomisticIGRefGCNSParameters(reference_activity=(3.0 / sf_V)u"Å^-3",
                                               species=:H)
        biased = MCAtomGrandCanonicalMoves(p_bias=0.4, bias_radius=2.0, bias_grid=8)
        @test_throws ArgumentError FreeBird.SamplingSchemes.nested_sampling_step!(ls2, params2, biased)
        galilean = MCAtomGrandCanonicalMoves(galilean_steps=2)
        @test_throws ArgumentError FreeBird.SamplingSchemes.nested_sampling_step!(ls2, params2, galilean)
    end

    function sf_run(seed; K=8, n_max=6, n_steps=60, tag="r", initialize=true, ls=nothing)
        Random.seed!(seed)
        ls = ls === nothing ? sf_liveset(K) : ls
        params = AtomisticIGRefGCNSParameters(mc_steps=40,
            reference_activity=(4.0 / sf_V)u"Å^-3", species=:H,
            allowed_fail_count=1000, n_max=n_max)
        df, lso, pout = ideal_gas_referenced_nested_sampling(
            ls, params, n_steps, MCAtomGrandCanonicalMoves(), sf_save(tag);
            initialize=initialize)
        sf_clean(tag)
        return df, lso, pout
    end

    @testset "surface driver end-to-end with same-seed determinism (seed 84004)" begin
        df1, lso1, _ = sf_run(84004)
        df2, lso2, _ = sf_run(84004)
        @test df1 == df2
        @test [ustrip(u"eV", w.energy) for w in lso1.walkers] ==
              [ustrip(u"eV", w.energy) for w in lso2.walkers]
        @test names(df1) == ["iter", "emax", "num_particles", "log_compression"]
        @test issorted(df1.emax, rev=true)
        @test maximum(df1.num_particles) <= 6
        @test maximum(w.list_num_par[1] for w in lso1.walkers) <= 6
        @test length(lso1.walkers) == 8
    end

    @testset "block continuation (initialize=false, seed 84005)" begin
        df1, lso1, _ = sf_run(84005; tag="c1")
        df2, lso2, _ = sf_run(84005; tag="c2", initialize=false, ls=lso1)
        @test nrow(df2) > 0
        @test df2.iter[1] == df1.iter[end] + 1
        @test df2.emax[1] <= df1.emax[end]
        @test issorted(df2.emax, rev=true)
        # a restored liveset violating the bounded cap is rejected
        @test any(w.list_num_par[1] > 2 for w in lso2.walkers)
        @test_throws ArgumentError sf_run(84005; tag="c3", initialize=false,
                                          ls=lso2, n_max=2)
        # the plain (non-surface) route continues the same way
        Random.seed!(84006)
        gen_ls = GenericAtomWalkers([AtomWalker{1}(sf_mkempty()) for _ in 1:8],
                                    LJParameters(epsilon=0.001, sigma=2.5, cutoff=1.8, shift=true))
        gp = AtomisticIGRefGCNSParameters(mc_steps=40,
            reference_activity=(4.0 / sf_V)u"Å^-3", species=:H, allowed_fail_count=1000)
        gdf1, glso, _ = ideal_gas_referenced_nested_sampling(
            gen_ls, gp, 40, MCAtomGrandCanonicalMoves(), sf_save("g1"))
        sf_clean("g1")
        gdf2, _, _ = ideal_gas_referenced_nested_sampling(
            glso, gp, 20, MCAtomGrandCanonicalMoves(), sf_save("g2"); initialize=false)
        sf_clean("g2")
        @test gdf2.iter[1] == gdf1.iter[end] + 1
        @test gdf2.emax[1] <= gdf1.emax[end]
    end

    @testset "capstone: bounded surface cross-route (seed 84621)" begin
        kb = 8.617333262e-5
        mass_ar = 39.948u"u"
        lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass_ar, T * u"K"))
        muS(zV, T) = (kb * T * (log(zV / sf_V) + 3 * log(lam(T)))) * u"eV"
        VqS = 1500.0u"Å^3"
        z0 = (4.0 / sf_V)u"Å^-3"
        n_cap = 4
        T0 = 300.0
        mus = [muS(2.0, T0), muS(4.0, T0)]

        Random.seed!(84621)
        ls = sf_liveset(16)
        params = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=z0,
                                              species=:H, allowed_fail_count=1000,
                                              n_max=n_cap)
        df, lso, _ = ideal_gas_referenced_nested_sampling(
            ls, params, 200, MCAtomGrandCanonicalMoves(), sf_save("cap"))
        sf_clean("cap")
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]
        @test maximum(df.num_particles) <= n_cap
        @test maximum(live_n) <= n_cap

        Random.seed!(85621)
        dfs = DataFrame[]
        lives = Vector{Float64}[]
        for N in 0:n_cap
            if N == 0
                push!(dfs, DataFrame(iter=Int[], emax=Float64[]))
                push!(lives, zeros(8))
                continue
            end
            walkers = AtomWalker{1}[]
            for _ in 1:16
                coords = [[rand() * 10.0, rand() * 10.0, rand() * 15.0] for _ in 1:N]
                push!(walkers, AtomWalker{1}(FastSystem(atomic_system(
                    [:H => SVector{3}(c)u"Å" for c in coords], sf_box, sf_pbc))))
            end
            lsN = LJSurfaceWalkers(walkers, sf_cps, sf_mksurf(); assign_energy=true)
            p = NestedSamplingParameters(mc_steps=40, initial_step_size=0.5, step_size=0.5,
                                         step_size_lo=0.01, step_size_up=2.0,
                                         allowed_fail_count=1000)
            dfN, lsoN, _ = nested_sampling(lsN, p, 120, MCRandomWalkClone(), sf_save("fx$N"))
            sf_clean("fx$N")
            push!(dfs, dfN)
            push!(lives, [ustrip(u"eV", w.energy) for w in lsoN.walkers])
        end

        ig = gc_thermodynamic_stats_ideal_ref(df, VqS, mass_ar, z0, mus, [T0]u"K";
                                              live_emax=live_e, live_numbers=live_n,
                                              n_max=n_cap)
        fx = gc_thermodynamic_stats_fixed_N(dfs, collect(0:n_cap), VqS, mass_ar,
                                            mus, [T0]u"K";
                                            n_walkers=16, live_emax=lives)
        @test maximum(abs.(ig.logXi .- fx.logXi)) < 1.0
        @test maximum(abs.(ig.mean_N .- fx.mean_N)) < 2.2
        @test minimum(ig.N_eff) > 3.5
    end
end
