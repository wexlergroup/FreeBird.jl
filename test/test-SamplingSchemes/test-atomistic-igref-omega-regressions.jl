# Cross-construction regression coverage for the chemical-potential-ordered
# atomistic ideal-gas-referenced route: an Omega-ordered run at its own mu,
# reduced over a temperature interval, against the energy-ordered run and the
# fixed-N stitch of canonical ladders (the independent oracle), plus the
# exactness identity at the residual-free temperature and the bitwise
# agreement of the effective-sample-size reduction on an Omega ledger.
#
# Calibration ledger (gates ship at >= 3x the maximum deviation over the seed
# pairs, per statistic; N_eff floors sit below the smallest observed value):
# - Surface three-way weld (the bounded surface fixture of the surface-route
#   testsets: z0V = 4, n_max = 4, K = 16; the Omega-ordered run at mu_run, the
#   chemical potential at which the ideal occupancy is 2 at 300 K, for 300
#   steps; the energy-ordered run for 200 steps; fixed-N ladders over N = 0..4
#   at K = 16 for 120 steps; targets (mu_run, T) for T in {250, 300, 400} K;
#   Omega-ordered igref seeds {1,2,3,86621}, energy-ordered igref seeds
#   {101,102,103,86721}, fixed-N seeds {1001,1002,1003,87621}). Every
#   Omega-ordered run ends on the empty-configuration atom (its live set all
#   empty), which at mu_run is the ground-state sector. Omega-ordered run
#   against the stitch: max devs logXi 0.672 (the shipped pair, 400 K) and
#   mean_N 0.417 (the seed-2 pair, 300 K; the shipped pair's largest is 0.329
#   at 300 K), N_eff floor observed 2.9 at 400 K; gates 2.1 and 1.3, floor 2.5.
#   Energy-ordered run against the stitch at 300 and 400 K: max devs logXi
#   0.480, mean_N 0.525 (400 K of the shipped pair, 300 K of the seed-2 pair),
#   N_eff floor observed 10.0; gates 1.5 and 1.6, floor 9. The energy-ordered
#   run's 250 K target lies outside its trust window (N_eff 1.2 to 4.3,
#   deviations up to 2.8 nats) and is deliberately not gated: at a fixed mu the
#   ideal occupancy falls to about 0.13 there, far from that run's reference
#   activity, whereas the Omega-ordered run, ordered at that mu, keeps N_eff
#   18 to 23 at 250 K. The residual-free temperature of this fixture's
#   reference activity is 0.1467 K (argon mass), where the exactness identity
#   below is evaluated.
# - Plain three-way weld (the 12 A periodic box of the ledger-assembly
#   testsets: z0V = 6, K = 24, 300 igref steps, fixed-N ladders over N = 0..12
#   at K = 24 for 150 steps; mu_run the chemical potential of ideal occupancy 4
#   at 300 K; Omega-ordered igref seeds {1,2,3,86631}, energy-ordered igref
#   seeds {101,102,103,86731}, fixed-N seeds {1001,1002,1003,87631}). At 300 K:
#   Omega-ordered run max devs logXi 0.211 (seed 3), mean_N 0.480 (seed 1),
#   N_eff floor observed 75.4; gates 0.65 and 1.5, floor 60. Energy-ordered run
#   max devs logXi 0.385 (the shipped seed), mean_N 0.792 (seed 1); gates 1.2
#   and 2.4. A 400 K target at this mu (ideal occupancy about 110, far above
#   the ladder's top at N = 12) lies outside the coverage of every construction
#   (deviations of 1 to 20 nats, N_eff 1 to 2) and is not part of the fixture.
# - The exactness identity and the effective-sample-size weld are floating-point
#   contracts and need no calibration.
@testset "cross-construction regressions for the chemical-potential-ordered route" begin
    using Random

    kb = 8.617333262e-5
    mass_ar = 39.948u"u"
    lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass_ar, T * u"K"))
    orsave(tag) = SaveEveryN(df_filename="_igoreg_$(tag).csv",
                             wk_filename="_igoreg_$(tag).traj.extxyz",
                             ls_filename="_igoreg_$(tag).ls.extxyz",
                             n_traj=10^7, n_snap=10^7, n_info=10^7)
    orclean(tag) = for f in ["_igoreg_$(tag).csv", "_igoreg_$(tag).traj.extxyz",
                             "_igoreg_$(tag).ls.extxyz"]
        rm(f, force=true)
    end

    # the surface fixture of the surface-route testsets
    sbox = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
    spbc = (true, true, false)
    sV = 1500.0
    sVq = 1500.0u"Å^3"
    surf_sys = FastSystem(atomic_system(
        [:H => [2.5, 2.5, 2.0]u"Å", :H => [7.5, 2.5, 2.0]u"Å",
         :H => [2.5, 7.5, 2.0]u"Å", :H => [7.5, 7.5, 2.0]u"Å"], sbox, spbc))
    mksurf() = AtomWalker(deepcopy(surf_sys); freeze_species=[:H])
    smkempty() = FastSystem(cell_vectors(surf_sys), periodicity(surf_sys),
                            empty(position(surf_sys, :)), empty(species(surf_sys, :)),
                            empty(mass(surf_sys, :)))
    cps = CompositeParameterSets(2, [LJParameters(epsilon=0.001, sigma=2.5, cutoff=1.8, shift=true),
                                     LJParameters(epsilon=0.003, sigma=2.5, cutoff=1.8, shift=true),
                                     LJParameters(epsilon=0.01, sigma=2.5, cutoff=1.8, shift=true)])
    sliveset(K) = LJSurfaceWalkers([AtomWalker{1}(smkempty()) for _ in 1:K], cps, mksurf();
                                   assign_energy=true)
    muS(zV, T) = (kb * T * (log(zV / sV) + 3 * log(lam(T)))) * u"eV"
    function surface_igref(seed, mu, n_steps; z0V=4.0, n_max=4, K=16)
        Random.seed!(seed)
        p = AtomisticIGRefGCNSParameters(mc_steps=40, reference_activity=(z0V / sV)u"Å^-3",
                                         species=:H, allowed_fail_count=1000, n_max=n_max,
                                         chemical_potential=mu)
        df, lso, _ = ideal_gas_referenced_nested_sampling(
            sliveset(K), p, n_steps, MCAtomGrandCanonicalMoves(), orsave("s"))
        orclean("s")
        return df, [ustrip(u"eV", w.energy) for w in lso.walkers], [w.list_num_par[1] for w in lso.walkers]
    end
    function surface_fixedN(seed; n_cap=4, K=16, n_steps=120)
        Random.seed!(seed)
        dfs = DataFrame[]
        lives = Vector{Float64}[]
        for N in 0:n_cap
            if N == 0
                # the N = 0 sector at zero energy: the four frozen H atoms sit 5 A
                # apart, beyond the 4.5 A cutoff, so the surface's own
                # energy_frozen_part is exactly zero and the empties of the igref
                # runs carry the same energy (a fixture with a nonzero frozen
                # part would shift this sector by e^{-beta E_frozen})
                push!(dfs, DataFrame(iter=Int[], emax=Float64[]))
                push!(lives, zeros(8))
                continue
            end
            walkers = AtomWalker{1}[]
            for _ in 1:K
                coords = [[rand() * 10.0, rand() * 10.0, rand() * 15.0] for _ in 1:N]
                push!(walkers, AtomWalker{1}(FastSystem(atomic_system(
                    [:H => SVector{3}(c)u"Å" for c in coords], sbox, spbc))))
            end
            lsN = LJSurfaceWalkers(walkers, cps, mksurf(); assign_energy=true)
            pN = NestedSamplingParameters(mc_steps=40, initial_step_size=0.5, step_size=0.5,
                                          step_size_lo=0.01, step_size_up=2.0,
                                          allowed_fail_count=1000)
            dfN, lsoN, _ = nested_sampling(lsN, pN, n_steps, MCRandomWalkClone(), orsave("fx$N"))
            orclean("fx$N")
            push!(dfs, dfN)
            push!(lives, [ustrip(u"eV", w.energy) for w in lsoN.walkers])
        end
        return dfs, lives
    end

    @testset "surface three-way weld over a temperature interval (seeds 86621, 86721, 87621)" begin
        z0s = (4.0 / sV)u"Å^-3"
        n_cap = 4
        Ts = [250.0, 300.0, 400.0]
        mu_run = muS(2.0, 300.0)
        df3, e3, n3 = surface_igref(86621, mu_run, 300)
        df2, e2, n2 = surface_igref(86721, 0.0u"eV", 200)
        dfs, lives = surface_fixedN(87621)
        # the Omega-ordered run ends on the empty-configuration atom, the
        # ground-state sector at this mu
        @test nrow(df3) > 0
        @test all(iszero, n3)
        @test issorted(df3.omega, rev=true)
        st3 = gc_thermodynamic_stats_ideal_ref(df3, sVq, mass_ar, z0s, [mu_run], Ts .* u"K";
                                               live_emax=e3, live_numbers=n3, n_max=n_cap)
        st2 = gc_thermodynamic_stats_ideal_ref(df2, sVq, mass_ar, z0s, [mu_run], Ts .* u"K";
                                               live_emax=e2, live_numbers=n2, n_max=n_cap)
        fx = gc_thermodynamic_stats_fixed_N(dfs, collect(0:n_cap), sVq, mass_ar, [mu_run], Ts .* u"K";
                                            n_walkers=16, live_emax=lives)
        # Omega-ordered run against the independent oracle at every target
        @test maximum(abs.(st3.logXi .- fx.logXi)) < 2.1
        @test maximum(abs.(st3.mean_N .- fx.mean_N)) < 1.3
        @test minimum(st3.N_eff) > 2.5
        # energy-ordered run against the oracle inside its own trust window
        @test maximum(abs.(st2.logXi[1, 2:3] .- fx.logXi[1, 2:3])) < 1.5
        @test maximum(abs.(st2.mean_N[1, 2:3] .- fx.mean_N[1, 2:3])) < 1.6
        @test minimum(st2.N_eff[1, 2:3]) > 9.0
        # at the residual-free temperature the reweighting factor at the run's
        # own mu is a function of Omega alone: a hand assembly on the :omega
        # column reproduces the reduction's logXi
        TL = reference_activity_temperature(z0s, mass_ar)
        stL = gc_thermodynamic_stats_ideal_ref(df3, sVq, mass_ar, z0s, [mu_run], [TL];
                                               live_emax=e3, live_numbers=n3, n_max=n_cap)
        beta = 1 / (kb * ustrip(u"K", TL))
        muv = ustrip(u"eV", mu_run)
        lc = Vector{Float64}(df3.log_compression)
        cs = cumsum(lc)
        log_w0 = vcat(0.0, cs[1:end-1]) .+ log.(-expm1.(lc))
        log_tail = cs[end] - log(length(e3))
        lw = vcat(log_w0 .- beta .* df3.omega,
                  fill(log_tail, length(e3)) .- beta .* (e3 .- muv .* n3))
        m = maximum(lw)
        hand = FreeBird.AnalysisTools._log_poisson_partial_sum(4.0, n_cap) + m + log(sum(exp.(lw .- m)))
        @test isapprox(hand, stL.logXi[1, 1]; rtol=1e-10)
        # the effective-sample-size reduction on an Omega ledger is the stats
        # N_eff, bitwise, and its relative form is finite over the interval
        ess = gc_effective_sample_size_ideal_ref(df3, sVq, mass_ar, z0s, [mu_run], Ts .* u"K";
                                                 live_emax=e3, live_numbers=n3)
        @test ess == st3.N_eff
        essr = gc_effective_sample_size_ideal_ref(df3, sVq, mass_ar, z0s, [mu_run], Ts .* u"K";
                                                  live_emax=e3, live_numbers=n3, relative=true)
        @test all(isfinite, essr)
    end

    @testset "plain three-way weld at 300 K (seeds 86631, 86731, 87631)" begin
        box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        pbc = (true, true, true)
        seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
        mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                               empty(position(seed_at, :)), empty(species(seed_at, :)),
                               empty(mass(seed_at, :)))
        V12 = 1728.0
        Vq12 = 1728.0u"Å^3"
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
        z0 = (6.0 / V12)u"Å^-3"
        mu12(zV, T) = (kb * T * (log(zV / V12) + 3 * log(lam(T)))) * u"eV"
        mu_run = mu12(4.0, 300.0)
        Ts = [300.0]
        function plain_igref(seed, mu, n_steps)
            Random.seed!(seed)
            p = AtomisticIGRefGCNSParameters(mc_steps=60, reference_activity=z0, species=:Ar,
                                             allowed_fail_count=100, chemical_potential=mu)
            df, lso, _ = ideal_gas_referenced_nested_sampling(
                GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:24], lj), p, n_steps,
                MCAtomGrandCanonicalMoves(), orsave("p"))
            orclean("p")
            return df, [ustrip(u"eV", w.energy) for w in lso.walkers], [w.list_num_par[1] for w in lso.walkers]
        end
        df3, e3, n3 = plain_igref(86631, mu_run, 300)
        df2, e2, n2 = plain_igref(86731, 0.0u"eV", 300)
        Random.seed!(87631)
        Nmax = 12
        dfs = DataFrame[]
        lives = Vector{Float64}[]
        for N in 0:Nmax
            if N == 0
                push!(dfs, DataFrame(iter=Int[], emax=Float64[]))
                push!(lives, zeros(8))
                continue
            end
            walkers = AtomWalker{1}[]
            for _ in 1:24
                coords = [[rand() * 12.0, rand() * 12.0, rand() * 12.0] for _ in 1:N]
                push!(walkers, AtomWalker{1}(FastSystem(atomic_system(
                    [:Ar => SVector{3}(c)u"Å" for c in coords], box, pbc))))
            end
            pN = NestedSamplingParameters(mc_steps=60, initial_step_size=0.5, step_size=0.5,
                                          step_size_lo=0.01, step_size_up=2.0,
                                          allowed_fail_count=1000)
            dfN, lsoN, _ = nested_sampling(GenericAtomWalkers(walkers, lj), pN, 150,
                                           MCRandomWalkClone(), orsave("px$N"))
            orclean("px$N")
            push!(dfs, dfN)
            push!(lives, [ustrip(u"eV", w.energy) for w in lsoN.walkers])
        end
        st3 = gc_thermodynamic_stats_ideal_ref(df3, Vq12, mass_ar, z0, [mu_run], Ts .* u"K";
                                               live_emax=e3, live_numbers=n3)
        st2 = gc_thermodynamic_stats_ideal_ref(df2, Vq12, mass_ar, z0, [mu_run], Ts .* u"K";
                                               live_emax=e2, live_numbers=n2)
        fx = gc_thermodynamic_stats_fixed_N(dfs, collect(0:Nmax), Vq12, mass_ar, [mu_run], Ts .* u"K";
                                            n_walkers=24, live_emax=lives)
        @test issorted(df3.omega, rev=true)
        @test abs(st3.logXi[1, 1] - fx.logXi[1, 1]) < 0.65
        @test abs(st3.mean_N[1, 1] - fx.mean_N[1, 1]) < 1.5
        @test st3.N_eff[1, 1] > 60.0
        @test abs(st2.logXi[1, 1] - fx.logXi[1, 1]) < 1.2
        @test abs(st2.mean_N[1, 1] - fx.mean_N[1, 1]) < 2.4
    end
end
