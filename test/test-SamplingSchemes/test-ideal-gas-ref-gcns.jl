@testset "Ideal-gas-referenced GC-NS" begin
    using Random

    kb = 8.617333262e-5  # eV/K

    # Shared 4x4 single-component square lattice template
    igref_L = 4
    igref_n_sites = igref_L * igref_L
    igref_template = MLattice{1,SquareLattice}(
        lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)],
        supercell_dimensions=(igref_L, igref_L, 1),
        periodicity=(true, true, false),
        cutoff_radii=[1.1, 1.5],
        components=[[false for _ in 1:igref_L*igref_L]],
        adsorptions=:full
    )

    function igref_run(template, ham, z0, n_walkers, n_steps; seed=7, mc_steps=100,
                       mc_routine=MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3))
        Random.seed!(seed)
        walkers = [LatticeWalker(deepcopy(template), energy=0.0u"eV", iter=0)
                   for _ in 1:n_walkers]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=mc_steps, reference_fugacity=z0)
        save = SaveEveryN("test_igref.csv", "test_igref.traj", "test_igref.ls",
                          100000, 100000, 100000)
        df, final_ls, params = ideal_gas_referenced_nested_sampling(
            liveset, params, Int64(n_steps), mc_routine, save)
        rm("test_igref.csv", force=true)
        rm("test_igref.traj", force=true)
        rm("test_igref.ls", force=true)
        live_E = [w.energy.val for w in final_ls.walkers]
        live_N = [sum(w.configuration.components[1]) for w in final_ls.walkers]
        return df, live_E, live_N
    end

    # ================================================================
    @testset "IdealGasReferencedGCNSParameters constructor" begin
        params = IdealGasReferencedGCNSParameters()
        @test params.mc_steps == 100
        @test params.reference_fugacity == 1.0
        @test params.energy_perturbation == 1e-12
        @test params.fail_count == 0
        @test params.allowed_fail_count == 10
        @test params.cluster_p == 0.3
        @test isempty(params.cluster_p_history)
        # move_stats defaults (issue #158)
        @test params.move_stats isa Dict{Symbol,Int}
        @test isempty(params.move_stats)

        params2 = IdealGasReferencedGCNSParameters(reference_fugacity=0.25, mc_steps=50)
        @test params2.reference_fugacity == 0.25
        @test params2.mc_steps == 50

        @test_throws ArgumentError IdealGasReferencedGCNSParameters(reference_fugacity=0.0)
        @test_throws ArgumentError IdealGasReferencedGCNSParameters(reference_fugacity=-1.0)
        # Tie-breaking is a correctness requirement on degenerate lattice spectra
        @test_throws ArgumentError IdealGasReferencedGCNSParameters(energy_perturbation=0.0)
    end

    # ================================================================
    @testset "nested_sampling_step! for IG-ref GC-NS" begin
        Random.seed!(3)
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        walkers = [LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
                   for _ in 1:5]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)

        params = IdealGasReferencedGCNSParameters(mc_steps=50, reference_fugacity=0.5)
        mc_routine = MCGrandCanonicalMoves()

        # Initialize walkers as draws from the Bernoulli(z0/(1+z0)) prior
        SamplingSchemes._init_ideal_gas_ref_walkers!(liveset, params)
        @test all(w.iter == 0 for w in liveset.walkers)

        e_type = typeof(walkers[1].energy)
        iter, emax, n_par, updated_liveset, updated_params = nested_sampling_step!(
            liveset, params, mc_routine)

        @test iter isa Union{Missing,Int}
        @test emax isa Union{Missing,e_type}
        @test n_par isa Union{Missing,Int}
        @test length(updated_liveset.walkers) == 5
        @test updated_params.fail_count >= 0
        if !(iter isa Missing)
            # The culled walker had the highest energy: all survivors sit
            # strictly below the recorded ceiling (perturbation breaks ties)
            @test all(w.energy < emax for w in updated_liveset.walkers)
        end
    end

    # ================================================================
    @testset "MC_grand_canonical_walk! z0 prior stationarity" begin
        # With a zero Hamiltonian (E ≡ 0) and a non-binding ceiling, the walk
        # must sample the Bernoulli(z0/(1+z0)) prior: ⟨N⟩ = M z0/(1+z0),
        # Var(N) = M z0/(1+z0)².
        ham0 = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        for z0 in (0.5, 2.0)
            Random.seed!(42)
            walker = LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
            MC_grand_canonical_walk!(500, walker, ham0, 1e3, 0.0;
                                     p_move=0.2, p_insert=0.4, z0=z0)
            acc_N = 0.0
            n_samples = 3000
            for _ in 1:n_samples
                MC_grand_canonical_walk!(10, walker, ham0, 1e3, 0.0;
                                         p_move=0.2, p_insert=0.4, z0=z0)
                acc_N += sum(walker.configuration.components[1])
            end
            p0 = z0 / (1 + z0)
            @test isapprox(acc_N / n_samples, igref_n_sites * p0; atol=0.5)
        end

        # ---- Biased-insertion stationarity grid (issue #158) ----
        # p_bias mixes uniform insertion with insertion restricted to
        # lattice_biased_sites(; predicate, shells=1); the composite MH
        # correction must leave the same Bernoulli(p0) product prior
        # invariant. p_bias = 0.5: the uniform component keeps the chain
        # irreducible, so a single-chain time average must reproduce
        # (N) = M p0, per-site occupancy p0, and Var(N) = M p0 (1 - p0).
        for (cell, (z0, pred)) in enumerate(Iterators.product((0.5, 2.0), (:contact, :cavity)))
            Random.seed!(4200 + cell)
            walker = LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
            MC_grand_canonical_walk!(300, walker, ham0, 1e3, 0.0;
                                     p_move=0.2, p_insert=0.4, z0=z0,
                                     p_bias=0.5, bias_predicate=pred, bias_shells=1)
            n_blocks = 1500
            Nsum = 0.0
            N2sum = 0.0
            site_sum = zeros(igref_n_sites)
            for _ in 1:n_blocks
                MC_grand_canonical_walk!(10, walker, ham0, 1e3, 0.0;
                                         p_move=0.2, p_insert=0.4, z0=z0,
                                         p_bias=0.5, bias_predicate=pred, bias_shells=1)
                occ = walker.configuration.components[1]
                Nsum += sum(occ)
                N2sum += sum(occ)^2
                site_sum .+= occ
            end
            p0 = z0 / (1 + z0)
            meanN = Nsum / n_blocks
            varN = N2sum / n_blocks - meanN^2
            # Tolerance model (10-step blocks, n = 1500 samples; Var(N) =
            # M p0 (1 - p0) = 3.56 at both z0): sigma((N)) ~ 0.10,
            # sigma(site) ~ 0.024 iid, SE(Var) ~ 0.26, all inflated by
            # residual block autocorrelation. Three-seed calibration (seed
            # bases 4200/5200/6200, all 4 cells): max |d(N)| = 0.145, max
            # site |d| = 0.052, max |dVar| = 0.361; shipped gates >= 3x the
            # maxima (0.5, 0.16, 1.2).
            @test isapprox(meanN, igref_n_sites * p0; atol=0.5)
            @test maximum(abs.(site_sum ./ n_blocks .- p0)) < 0.16
            @test isapprox(varN, igref_n_sites * p0 * (1 - p0); atol=1.2)
        end

        # p_bias = 1.0: the kernel is still prior-invariant (each sub-kernel
        # is in detailed balance) but REDUCIBLE: :contact cannot leave N = 0
        # and :cavity can never insert next to the occupied region, so dense
        # states are unreachable and time averages are meaningless (hence the
        # constructor warning and the non-ergodicity pin below). Invariance is
        # what the composite MH factor must guarantee, so test it directly: an
        # ensemble of chains started from EXACT prior draws must remain
        # prior-distributed after any number of steps.
        for (cell, (z0, pred)) in enumerate(Iterators.product((0.5, 2.0), (:contact, :cavity)))
            Random.seed!(4300 + cell)
            p0 = z0 / (1 + z0)
            n_chains = 1200
            Nsum = 0.0
            N2sum = 0.0
            site_sum = zeros(igref_n_sites)
            for _ in 1:n_chains
                walker = LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
                random_microstate!(walker.configuration; p=p0)   # exact prior draw
                MC_grand_canonical_walk!(25, walker, ham0, 1e3, 0.0;
                                         p_move=0.2, p_insert=0.4, z0=z0,
                                         p_bias=1.0, bias_predicate=pred, bias_shells=1)
                occ = walker.configuration.components[1]
                Nsum += sum(occ)
                N2sum += sum(occ)^2
                site_sum .+= occ
            end
            meanN = Nsum / n_chains
            varN = N2sum / n_chains - meanN^2
            # Chains are iid (fresh prior starts): sigma((N)) = 0.054,
            # sigma(site) = 0.0136, SE(Var) ~ 0.145 at n = 1200. Three-seed
            # calibration (seed bases 4300/5300/6300, all 4 cells): max
            # |d(N)| = 0.083, max site |d| = 0.038, max |dVar| = 0.272;
            # shipped gates >= 3x the maxima (0.35, 0.12, 0.85).
            @test isapprox(meanN, igref_n_sites * p0; atol=0.35)
            @test maximum(abs.(site_sum ./ n_chains .- p0)) < 0.12
            @test isapprox(varN, igref_n_sites * p0 * (1 - p0); atol=0.85)
        end

        # z0 must be positive
        walker = LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
        @test_throws ArgumentError MC_grand_canonical_walk!(
            1, walker, ham0, 1e3, 0.0; z0=0.0)
    end

    # ================================================================
    @testset "Exact microstate distribution on 2x2 (biased kernel)" begin
        # 2x2 torus: each site's 4 NN collapse onto 2 distinct sites, so the
        # default min-image construction warns; multiplicity lists are exact
        # here and the bias predicates only read occupancy, so duplicated
        # neighbor entries are harmless.
        tiny_template = MLattice{1,SquareLattice}(
            lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(2, 2, 1), periodicity=(true, true, false),
            cutoff_radii=[1.1], components=[[false for _ in 1:4]],
            adsorptions=:full, image_multiplicity=true)
        ham0 = GenericLatticeHamiltonian(0.0, [0.0], u"eV")
        z0 = 2.0
        # Zero Hamiltonian, non-binding ceiling: the walk must sample the
        # Bernoulli product prior exactly, P(state) = z0^N / (1+z0)^4
        for (ip, pred) in enumerate((:contact, :cavity))
            Random.seed!(5150 + ip)
            walker = LatticeWalker(deepcopy(tiny_template), energy=0.0u"eV", iter=0)
            MC_grand_canonical_walk!(200, walker, ham0, 1e3, 0.0;
                                     p_move=0.2, p_insert=0.4, z0=z0,
                                     p_bias=0.5, bias_predicate=pred, bias_shells=1)
            n_samples = 20_000
            counts = zeros(Int, 16)
            Nsum = 0.0
            for _ in 1:n_samples
                # 8-step blocks (~2 sweeps of the 4-site lattice) decorrelate
                MC_grand_canonical_walk!(8, walker, ham0, 1e3, 0.0;
                                         p_move=0.2, p_insert=0.4, z0=z0,
                                         p_bias=0.5, bias_predicate=pred, bias_shells=1)
                occ = walker.configuration.components[1]
                mask = sum(Int(occ[s]) << (s - 1) for s in 1:4)
                counts[mask + 1] += 1
                Nsum += sum(occ)
            end
            # Multinomial gates: per cell sigma = sqrt(p(1-p)/n) at n = 20000;
            # the max over 16 cells x 2 predicates of ~3 sigma is the
            # expected extreme of 32 comparisons, inflated by residual block
            # autocorrelation. Three-seed calibration (seed bases 5150/6150/
            # 7150): max |dev| = 3.0 iid-sigma, max |d(N)| = 0.023; shipped
            # gates >= 3x the maxima (k = 9, atol 0.07). Discrimination: a
            # deliberately broken kernel (reverse density on the pre-delete
            # configuration, or the uniform-proposal ratio on biased draws)
            # fails these gates at 20-55 sigma on the discriminating cells;
            # per-state residuals below ~5e-3 sit under them, the honest
            # floor at feasible n.
            for mask in 0:15
                p = z0^count_ones(mask) / (1 + z0)^4
                sig = sqrt(p * (1 - p) / n_samples)
                @test abs(counts[mask + 1] / n_samples - p) < 9 * sig
            end
            @test isapprox(Nsum / n_samples, 4 * z0 / (1 + z0); atol=0.07)
        end
    end

    # ================================================================
    @testset "p_bias = 1.0 :contact is non-ergodic from empty" begin
        # From N = 0 the contact set is empty, the biased branch has no
        # proposal, and the uniform branch is never selected: N == 0 is an
        # exact invariant even at z0 = 2 insertion pressure. This is the
        # reducibility the MCGrandCanonicalMoves constructor warns about.
        ham0_ne = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        Random.seed!(158)
        walker = LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
        for _ in 1:10
            MC_grand_canonical_walk!(200, walker, ham0_ne, 1e3, 0.0;
                p_move=0.2, p_insert=0.4, z0=2.0,
                p_bias=1.0, bias_predicate=:contact, bias_shells=1)
            @test sum(walker.configuration.components[1]) == 0
        end
    end

    # ================================================================
    @testset "insert branch split matches p_bias (counters)" begin
        ham0_bs = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        Random.seed!(31)
        lat_bs = deepcopy(igref_template)
        lat_bs.components[1][1:8] .= true
        walker = LatticeWalker(lat_bs, energy=0.0u"eV", iter=0)
        _, _, _, _, _, c = MC_grand_canonical_walk!(
            6000, walker, ham0_bs, 1e3, 0.0;
            p_move=0.2, p_insert=0.4, z0=1.0,
            p_bias=0.5, bias_predicate=:contact, bias_shells=1)
        k = c.insert_biased_attempted
        n = k + c.insert_uniform_attempted
        # n ~ 0.4 x 6000 = 2400 insert-branch trials, branch coin
        # Bernoulli(1/2): k ~ Binomial(n, 1/2), sigma = sqrt(n/4) ~ 24.5;
        # gate 4 sigma (~98 counts, alpha ~ 6e-5). :contact skips only at
        # N = 0 (prior weight 2^-16 at z0 = 1: well under one expected
        # skipped trial over the run), negligible against the gate.
        @test n > 2000
        @test abs(k - n / 2) <= 4 * sqrt(n / 4)
    end

    # ================================================================
    @testset "Langmuir closed form (non-interacting reference)" begin
        # Zero neighbor interactions, on-site adsorption only: E = ε_ads·N, so
        # Ξ(μ,T) = (1 + z·e^{-βε})^M exactly (Langmuir isotherm) — the lattice
        # analog of the atomistic IdealGasParameters validation.
        eps_ads = -0.04
        ham_ni = GenericLatticeHamiltonian(eps_ads, [0.0, 0.0], u"eV")
        n_walkers = 100

        df, live_E, live_N = igref_run(igref_template, ham_ni, 1.0, n_walkers, 3000)
        @test nrow(df) > 0
        @test names(df) == ["iter", "emax", "num_particles"]
        # E-sorted NS: recorded energy ceilings are non-increasing
        @test issorted(df.emax, rev=true)

        # μ grid spans the reliable reweighting window around μ_ref(T) = kT·ln z0 = 0
        # plus one deliberately out-of-window point (μ = -0.08) for the N_eff check
        μs = [-0.08, -0.04, -0.03, -0.02]
        Ts = [200.0, 300.0, 500.0]
        stats = gc_thermodynamic_stats_ideal_ref(
            df, igref_n_sites, 1.0, μs, Ts, n_walkers;
            ω0=(n_walkers + 1) / n_walkers, live_emax=live_E, live_numbers=live_N)
        stats_notail = gc_thermodynamic_stats_ideal_ref(
            df, igref_n_sites, 1.0, μs, Ts, n_walkers)

        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            i == 1 && continue  # skip the out-of-window point
            β = 1 / (kb * T)
            x = exp(β * μ) * exp(-β * eps_ads)
            logXi_exact = igref_n_sites * log1p(x)
            meanN_exact = igref_n_sites * x / (1 + x)
            @test isapprox(stats.logXi[i, j], logXi_exact; atol=0.5)
            @test isapprox(stats.mean_N[i, j], meanN_exact; rtol=0.10)
            @test isapprox(stats.mean_U[i, j], eps_ads * meanN_exact; rtol=0.10)
        end

        # The live-set tail adds positive prior mass: logΞ strictly increases
        @test all(stats.logXi .> stats_notail.logXi)

        # N_eff collapses out of the reweighting window: the far point
        # (μ = -0.08 at T = 200 K, |βμ - ln z0| ≈ 4.6) must have a much
        # smaller effective sample size than an in-window point
        @test stats.N_eff[1, 1] < stats.N_eff[4, 1] / 10
        @test all(stats.N_eff .<= nrow(df) + length(live_E))

        # The effective-sample-size reduction reproduces the collapse
        # ordering and the row-count ceiling on the same fixture
        ess_l = gc_effective_sample_size_ideal_ref(
            df, igref_n_sites, 1.0, μs, Ts, n_walkers;
            ω0=(n_walkers + 1) / n_walkers, live_emax=live_E, live_numbers=live_N)
        @test ess_l[1, 1] < ess_l[4, 1] / 10
        @test all(ess_l .<= nrow(df) + length(live_E))
    end

    # ================================================================
    @testset "Interacting lattice vs exact enumeration (z0 ≠ 1)" begin
        # NN + NNN attractive lattice gas; reference fugacity centered on the
        # target (μ = -0.05 eV, T = 300 K), i.e. z0 = exp(βμ) ≈ 0.145. This
        # exercises the z0-weighted insert/delete ratios, the Bernoulli(p0)
        # initialization, and the (1+z0)^M absolute normalization end to end.
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        μ_center = -0.05
        T_center = 300.0
        z0 = exp(μ_center / (kb * T_center))
        n_walkers = 100

        df, live_E, live_N = igref_run(igref_template, ham, z0, n_walkers, 4000; seed=11)
        @test nrow(df) > 0
        @test issorted(df.emax, rev=true)
        # The energy ladder must reach deep order (ground state E = -1.04 eV:
        # full lattice, 16·(-0.04) + 32·(-0.01) + 32·(-0.0025))
        @test minimum(df.emax) < -1.0

        μs = [-0.06, -0.05]
        Ts = [300.0, 400.0]

        # Exact grand-canonical reference: enumerate all 2^16 microstates once
        E_vals = Vector{Float64}(undef, 2^igref_n_sites)
        N_vals = Vector{Int}(undef, 2^igref_n_sites)
        lattice = deepcopy(igref_template)
        for mask in 0:(2^igref_n_sites - 1)
            for site in 1:igref_n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_vals[mask+1] = interacting_energy(lattice, ham).val
            N_vals[mask+1] = sum(lattice.components[1])
        end

        stats = gc_thermodynamic_stats_ideal_ref(
            df, igref_n_sites, z0, μs, Ts, n_walkers;
            ω0=(n_walkers + 1) / n_walkers, live_emax=live_E, live_numbers=live_N)

        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            β = 1 / (kb * T)
            log_w = β .* μ .* N_vals .- β .* E_vals
            max_log = maximum(log_w)
            w = exp.(log_w .- max_log)
            Z = sum(w)
            logXi_exact = max_log + log(Z)
            meanN_exact = sum(w .* N_vals) / Z
            meanU_exact = sum(w .* E_vals) / Z
            varN_exact = sum(w .* N_vals .^ 2) / Z - meanN_exact^2

            # Tolerances: logΞ atol = 1.5 ≈ 2.7x the observed worst
            # seed-to-seed error (±0.55 over 6 seeds), kept below the ≈2.2
            # shift a missing (1+z0)^M normalization would cause so the test
            # still discriminates; ⟨N⟩/⟨U⟩ rtol = 0.10 ≈ 4x the ±0.27 spread.
            # This is the first absolute-Ξ check in the test suite.
            @test isapprox(stats.logXi[i, j], logXi_exact; atol=1.5)
            @test isapprox(stats.mean_N[i, j], meanN_exact; rtol=0.10)
            @test isapprox(stats.mean_U[i, j], meanU_exact; rtol=0.10)
            @test isapprox(stats.var_N[i, j], varN_exact; rtol=0.25)
        end
    end

    # ================================================================
    @testset "Interacting lattice vs exact enumeration (biased :contact)" begin
        # Same physics as the unbiased z0 != 1 testset above, driven through
        # the composite kernel with p_bias = 0.5 :contact via the NS step:
        # the end-to-end check that biased insertion leaves the ideal-gas
        # reference bookkeeping (weights, (1+z0)^M normalization) intact.
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        μ_center = -0.05
        T_center = 300.0
        z0 = exp(μ_center / (kb * T_center))
        n_walkers = 100

        mc_bias = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, p_bias=0.5)
        df, live_E, live_N = igref_run(igref_template, ham, z0, n_walkers, 4000;
                                       seed=13, mc_routine=mc_bias)
        @test nrow(df) > 0
        @test issorted(df.emax, rev=true)
        # The ladder reaches deep order (ground state -1.04 eV)
        @test minimum(df.emax) < -1.0

        # mu grid: two in-window points plus one deliberately out-of-window
        # point (mu = -0.08) for the Kish N_eff gate
        μs = [-0.08, -0.06, -0.05]
        Ts = [300.0, 400.0]

        E_vals = Vector{Float64}(undef, 2^igref_n_sites)
        N_vals = Vector{Int}(undef, 2^igref_n_sites)
        lattice = deepcopy(igref_template)
        for mask in 0:(2^igref_n_sites - 1)
            for site in 1:igref_n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_vals[mask+1] = interacting_energy(lattice, ham).val
            N_vals[mask+1] = sum(lattice.components[1])
        end

        stats = gc_thermodynamic_stats_ideal_ref(
            df, igref_n_sites, z0, μs, Ts, n_walkers;
            ω0=(n_walkers + 1) / n_walkers, live_emax=live_E, live_numbers=live_N)

        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            i == 1 && continue  # the out-of-window point feeds only the N_eff gate
            β = 1 / (kb * T)
            log_w = β .* μ .* N_vals .- β .* E_vals
            max_log = maximum(log_w)
            w = exp.(log_w .- max_log)
            Z = sum(w)
            logXi_exact = max_log + log(Z)
            meanN_exact = sum(w .* N_vals) / Z
            meanU_exact = sum(w .* E_vals) / Z
            varN_exact = sum(w .* N_vals .^ 2) / Z - meanN_exact^2
            # Tolerances recalibrated for the biased kernel over 3 seeds
            # (13/14/15; every in-window point gated at N_eff >= 1747): max
            # |d logXi| = 0.45, max relative d(N) = 0.044, d(U) = 0.052;
            # shipped >= 3x the maxima (1.5, 0.15, 0.16), still below the
            # ~2.2 shift a missing (1+z0)^M normalization would cause. Var N
            # is not gated here: at these near-step grid points its
            # reweighting noise reaches 0.28 relative (a >= 3x gate would be
            # vacuous), and it is exercised by the unbiased register above.
            @test isapprox(stats.logXi[i, j], logXi_exact; atol=1.5)
            @test isapprox(stats.mean_N[i, j], meanN_exact; rtol=0.15)
            @test isapprox(stats.mean_U[i, j], meanU_exact; rtol=0.16)
        end

        # Kish N_eff degrades out of the reweighting window; every reweighted
        # point stays bounded by the sample count
        @test stats.N_eff[1, 1] < stats.N_eff[3, 1]
        @test all(stats.N_eff .<= nrow(df) + length(live_E))
    end

    # ================================================================
    @testset "Exact Skilling closure (E ≡ 0)" begin
        # With a zero Hamiltonian every configuration has E = 0 (up to the
        # 1e-12 tie-breaking tags), so at z = z0 the reweighting factor is 1
        # and logΞ must equal M·ln(1+z0) EXACTLY when the dead weights use
        # ω0 = (K+n_cull)/K and the live-set tail closes the ladder:
        # Σω = (1 - r^n) + r^n = 1, independent of n and of the sampled N_j.
        # This is an algebraic identity test of the weight bookkeeping — it
        # would catch any ω0/tail double counting.
        ham0 = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        n_walkers = 50
        df, live_E, live_N = igref_run(igref_template, ham0, 1.0, n_walkers, 300;
                                       seed=17, mc_steps=30)
        @test nrow(df) > 0
        stats = gc_thermodynamic_stats_ideal_ref(
            df, igref_n_sites, 1.0, [0.0], [300.0], n_walkers;
            ω0=(n_walkers + 1) / n_walkers, live_emax=live_E, live_numbers=live_N)
        # β|E| ≲ 4e-11 from the tie-breaking tags; atol dominated by that
        @test isapprox(stats.logXi[1, 1], igref_n_sites * log(2.0); atol=1e-6)
    end

    # ================================================================
    @testset "gc_thermodynamic_stats_ideal_ref argument validation" begin
        df_ok = DataFrame(iter=[1, 2], emax=[-0.1, -0.2], num_particles=[3, 4])
        μs = [-0.05]
        Ts = [300.0]

        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_ok, 16, 0.0, μs, Ts, 100)
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_ok, 0, 1.0, μs, Ts, 100)
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_ok, 16, 1.0, μs, Ts, 100; live_emax=[-0.1])
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_ok, 16, 1.0, μs, Ts, 100; live_numbers=[1])
        @test_throws DimensionMismatch gc_thermodynamic_stats_ideal_ref(
            df_ok, 16, 1.0, μs, Ts, 100; live_emax=[-0.1], live_numbers=[1, 2])

        df_empty = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[])
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_empty, 16, 1.0, μs, Ts, 100)

        # Live walkers alone are a valid input (a run that never culled:
        # n_iters = 0, each live walker carries weight 1/K) — check the
        # closed form logΞ = M·ln(1+z0) + logsumexp(-ln K + N·βμ - β·E)
        live_E = [-0.1, -0.2]
        live_N = [1, 2]
        stats_live = @test_logs (:warn, r"live-set tail") match_mode = :any gc_thermodynamic_stats_ideal_ref(
            df_empty, 16, 1.0, μs, Ts, 100;
            live_emax=live_E, live_numbers=live_N)
        β = 1 / (kb * Ts[1])
        terms = [exp(-log(100) + live_N[k] * β * μs[1] - β * live_E[k]) for k in 1:2]
        logXi_expected = 16 * log(2.0) + log(sum(terms))
        @test isapprox(stats_live.logXi[1, 1], logXi_expected; rtol=1e-10)
        @test isapprox(stats_live.mean_N[1, 1],
                       sum(terms .* live_N) / sum(terms); rtol=1e-10)

        # Output shape
        stats = gc_thermodynamic_stats_ideal_ref(
            df_ok, 16, 1.0, [-0.05, -0.04, -0.03], [200.0, 300.0], 100)
        @test size(stats.logXi) == (3, 2)
        @test size(stats.var_N) == (3, 2)
        @test size(stats.N_eff) == (3, 2)
    end

    # ================================================================
    @testset "Cluster moves through the IG-ref step" begin
        # Exercise the clusters_freq > 0 path: geometric cluster moves inside
        # MC_grand_canonical_walk! plus the adaptive cluster_p machinery on
        # IdealGasReferencedGCNSParameters
        Random.seed!(21)
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        walkers = [LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
                   for _ in 1:20]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=50, reference_fugacity=0.5)
        mc_routine = MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25,
                                           clusters_freq=1, swaps_freq=1)
        save = SaveEveryN("test_igref_cl.csv", "test_igref_cl.traj", "test_igref_cl.ls",
                          100000, 100000, 100000)
        df, _, params_out = ideal_gas_referenced_nested_sampling(
            liveset, params, Int64(400), mc_routine, save)
        rm("test_igref_cl.csv", force=true)
        rm("test_igref_cl.traj", force=true)
        rm("test_igref_cl.ls", force=true)

        @test nrow(df) > 0
        # Adaptation must have fired: history populated, cluster_p in bounds
        @test !isempty(params_out.cluster_p_history)
        @test all(mc_routine.cluster_p_floor .<= params_out.cluster_p_history
                  .<= mc_routine.cluster_p_ceiling)
        @test length(params_out.cluster_accept_history) ==
              length(params_out.cluster_p_history)
    end

    # ================================================================
    @testset "Same-seed reproducibility" begin
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        df1, live_E1, _ = igref_run(igref_template, ham, 2.0, 30, 300; seed=99, mc_steps=50)
        df2, live_E2, _ = igref_run(igref_template, ham, 2.0, 30, 300; seed=99, mc_steps=50)
        @test df1.emax == df2.emax
        @test df1.num_particles == df2.num_particles
        @test live_E1 == live_E2
    end

    # ================================================================
    @testset "Same-seed A/B: p_bias = 0.0 is RNG-identical to default" begin
        # The p_bias = 0 code path must draw exactly the same RNG stream as
        # a routine that never mentions p_bias; any extra rand() would
        # silently break reproducibility of published runs.
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        mcA = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)
        mcB = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, p_bias=0.0)
        dfA, live_EA, _ = igref_run(igref_template, ham, 2.0, 30, 300;
                                    seed=99, mc_steps=50, mc_routine=mcA)
        dfB, live_EB, _ = igref_run(igref_template, ham, 2.0, 30, 300;
                                    seed=99, mc_steps=50, mc_routine=mcB)
        @test dfA.emax == dfB.emax
        @test dfA.num_particles == dfB.num_particles
        @test live_EA == live_EB
    end

    # ================================================================
    @testset "dead-point callback and energy substitution (IG-ref route)" begin
        using Random
        dpc_lat = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:16]],
            adsorptions=:full)
        dpc_ham_a = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        dpc_ham_b = GenericLatticeHamiltonian(-0.04, [-0.02, -0.0025], u"eV")
        dpc_save = SaveEveryN("t_dpc_ig.csv", "t_dpc_ig.traj", "t_dpc_ig.ls",
                              1000000, 1000000, 1000000)
        dpc_cleanup() = rm.(["t_dpc_ig.csv", "t_dpc_ig.traj", "t_dpc_ig.ls"],
                            force=true)

        # Exact grand sums for both Hamiltonians by per-sector enumeration
        # over all 2^16 configurations, computed before any seeding because
        # exact_enumeration consumes the global RNG
        dpc_kb = 8.617333262e-5
        dpc_T = 600.0
        dpc_mu = 0.0
        dpc_beta = 1 / (dpc_kb * dpc_T)
        function dpc_exact_lnXi(h)
            lnZs = Float64[]
            for N in 0:16
                occ = vcat(fill(true, N), fill(false, 16 - N))
                latN = deepcopy(dpc_lat)
                latN.components[1] .= occ
                dfe, _ = exact_enumeration(latN, h)
                Es = [ustrip(u"eV", e) for e in dfe.energy]
                Emin = minimum(Es)
                push!(lnZs, dpc_beta * dpc_mu * N - dpc_beta * Emin +
                            log(sum(exp.(-dpc_beta .* (Es .- Emin)))))
            end
            m = maximum(lnZs)
            return m + log(sum(exp.(lnZs .- m)))
        end
        lnXi_a = dpc_exact_lnXi(dpc_ham_a)
        lnXi_b = dpc_exact_lnXi(dpc_ham_b)

        # One seeded ladder under Hamiltonian A, configurations collected
        # through the callback; K = 64, full descent against the
        # 16·ln 2 ≈ 11.1 nat reference depth
        dpc_K = 64
        Random.seed!(4371)
        walkers = [LatticeWalker(deepcopy(dpc_lat), energy=0.0u"eV", iter=0)
                   for _ in 1:dpc_K]
        ls = LatticeGasWalkers(walkers, dpc_ham_a; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=40,
            reference_fugacity=1.0, energy_perturbation=1e-9)
        configs = Vector{Bool}[]
        seen = Tuple{Int,Float64,Int}[]
        df, ls_out, _ = ideal_gas_referenced_nested_sampling(ls, params,
            Int64(820), MCGrandCanonicalMoves(), dpc_save;
            dead_point_callback=(iter, w) -> begin
                push!(configs, copy(w.configuration.components[1]))
                push!(seen, (iter, w.energy.val,
                             sum(w.configuration.components[1])))
            end)
        dpc_cleanup()

        # Invocation count equals nrow(df); the (iter, energy, N) triple seen
        # by the callback matches the ledger row bit-exactly
        @test nrow(df) > 0
        @test length(seen) == nrow(df)
        @test [t[1] for t in seen] == df.iter
        @test [t[2] for t in seen] == df.emax
        @test [t[3] for t in seen] == df.num_particles

        # Energy faithfulness: recomputing Hamiltonian A on a collected
        # configuration matches the ledger energy within the perturbation
        # half-width (uniform on ±energy_perturbation/2)
        dpc_probe = deepcopy(dpc_lat)
        function dpc_recompute(h, c)
            dpc_probe.components[1] .= c
            return ustrip(u"eV", interacting_energy(dpc_probe, h))
        end
        @test all(abs(dpc_recompute(dpc_ham_a, c) - e) <= 5.1e-10
                  for (c, e) in zip(configs, df.emax))

        # Substitution closure: Hamiltonian B's energies, recomputed offline
        # from the collected configurations, replace the ledger column and
        # the live tail, and gc_thermodynamic_stats_ideal_ref then closes
        # against Hamiltonian B's exact grand sum with the SAME Skilling
        # weights; the control route closes against Hamiltonian A's.
        # Tolerances calibrated at seeds 4371/4372/4373: max |dev| 0.236 (A)
        # and 0.343 (B) against σ = √(16·ln 2/64) ≈ 0.42; shipped at ≥ 3×.
        live_Ea = [w.energy.val for w in ls_out.walkers]
        live_N = [Int(sum(w.configuration.components[1])) for w in ls_out.walkers]
        sa = gc_thermodynamic_stats_ideal_ref(df, 16, 1.0, [dpc_mu], [dpc_T],
            dpc_K; ω0=(dpc_K + 1) / dpc_K, live_emax=live_Ea, live_numbers=live_N)
        df_b = copy(df)
        df_b.emax = [dpc_recompute(dpc_ham_b, c) for c in configs]
        live_Eb = [dpc_recompute(dpc_ham_b, w.configuration.components[1])
                   for w in ls_out.walkers]
        sb = gc_thermodynamic_stats_ideal_ref(df_b, 16, 1.0, [dpc_mu], [dpc_T],
            dpc_K; ω0=(dpc_K + 1) / dpc_K, live_emax=live_Eb, live_numbers=live_N)
        @test abs(sa.logXi[1, 1] - lnXi_a) < 0.75
        @test abs(sb.logXi[1, 1] - lnXi_b) < 1.1
        # The substituted estimator's Kish N_eff is a real diagnostic: finite,
        # positive, and no larger than the control's at the same grid point
        @test 0 < sb.N_eff[1, 1] <= sa.N_eff[1, 1]

        # Stream neutrality: same-seed A/B with and without the callback
        function dpc_ab(seed, cb)
            Random.seed!(seed)
            ws = [LatticeWalker(deepcopy(dpc_lat), energy=0.0u"eV", iter=0)
                  for _ in 1:16]
            l = LatticeGasWalkers(ws, dpc_ham_a; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=20,
                reference_fugacity=1.0, energy_perturbation=1e-9)
            d, _, _ = ideal_gas_referenced_nested_sampling(l, p, Int64(100),
                MCGrandCanonicalMoves(), dpc_save; dead_point_callback=cb)
            dpc_cleanup()
            return d
        end
        dfA = dpc_ab(4373, nothing)
        dfB = dpc_ab(4373, (iter, w) -> nothing)
        @test dfA.iter == dfB.iter
        @test dfA.emax == dfB.emax
        @test dfA.num_particles == dfB.num_particles
    end

    @testset "opt-in incremental walk" begin
        # Calibration ledger (gates at >= 3x the max three-seed deviation,
        # per statistic): incremental = true full-descent ladders on the
        # 4x4 two-shell fixture (K = 64, mc_steps = 40, 820 steps, z0 = 1,
        # reduced at mu = 0, T = 600 K), seeds 96001/96002/96003 measured
        # |d logXi| = 0.153 / 0.282 / 0.156 (sigma class sqrt(16 ln2 / 64)
        # ~ 0.42, matching the shipped substitution-closure calibration) and
        # |d mean_N| = 0.057 / 0.088 / 0.076: gates 0.85 and 0.27.
        inc_lat() = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:16]], adsorptions=:full)
        inc_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        inc_save = SaveEveryN("t_inc.csv", "t_inc.traj", "t_inc.ls",
                              1000000, 1000000, 1000000)
        inc_cleanup() = rm.(["t_inc.csv", "t_inc.traj", "t_inc.ls"],
                            force=true)
        inc_T = 600.0
        inc_mu = 0.0
        inc_beta = 1 / (kb * inc_T)

        # Exact grand sum by per-sector enumeration, computed before any
        # seeding because exact_enumeration consumes the global RNG
        inc_lnZs = Float64[]
        for N in 0:16
            occ = vcat(fill(true, N), fill(false, 16 - N))
            latN = inc_lat()
            latN.components[1] .= occ
            dfe, _ = exact_enumeration(latN, inc_ham)
            Es = [ustrip(u"eV", e) for e in dfe.energy]
            Emin = minimum(Es)
            push!(inc_lnZs, inc_beta * inc_mu * N - inc_beta * Emin +
                            log(sum(exp.(-inc_beta .* (Es .- Emin)))))
        end
        inc_m = maximum(inc_lnZs)
        inc_lnXi = inc_m + log(sum(exp.(inc_lnZs .- inc_m)))
        inc_meanN = sum((0:16) .* exp.(inc_lnZs .- inc_m)) /
                    sum(exp.(inc_lnZs .- inc_m))

        # Same-seed A/B stream identity: a routine that never mentions the
        # field against incremental = false, on the driver
        function inc_ab(seed, routine)
            Random.seed!(seed)
            ws = [LatticeWalker(deepcopy(inc_lat()), energy=0.0u"eV", iter=0)
                  for _ in 1:16]
            ls = LatticeGasWalkers(ws, inc_ham; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=20,
                reference_fugacity=1.0, energy_perturbation=1e-9)
            d, lsx, _ = ideal_gas_referenced_nested_sampling(ls, p,
                Int64(100), routine, inc_save)
            inc_cleanup()
            return d, lsx
        end
        dA, lsA = inc_ab(97001, MCGrandCanonicalMoves())
        dB, lsB = inc_ab(97001, MCGrandCanonicalMoves(incremental=false))
        @test dA.iter == dB.iter
        @test dA.emax == dB.emax
        @test dA.num_particles == dB.num_particles
        @test [w.energy.val for w in lsA.walkers] ==
              [w.energy.val for w in lsB.walkers]

        # Statistical gate against the exact grand sum, three seeds
        for seed in (96001, 96002, 96003)
            Random.seed!(seed)
            ws = [LatticeWalker(deepcopy(inc_lat()), energy=0.0u"eV", iter=0)
                  for _ in 1:64]
            ls = LatticeGasWalkers(ws, inc_ham; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=40,
                reference_fugacity=1.0, energy_perturbation=1e-9)
            d, lsx, _ = ideal_gas_referenced_nested_sampling(ls, p,
                Int64(820), MCGrandCanonicalMoves(incremental=true), inc_save)
            inc_cleanup()
            live_E = [w.energy.val for w in lsx.walkers]
            live_N = [Int(sum(w.configuration.components[1]))
                      for w in lsx.walkers]
            s = gc_thermodynamic_stats_ideal_ref(d, 16, 1.0, [inc_mu],
                [inc_T], 64; ω0=65 / 64, live_emax=live_E,
                live_numbers=live_N)
            @test abs(s.logXi[1, 1] - inc_lnXi) < 0.85
            @test abs(s.mean_N[1, 1] - inc_meanN) < 0.27
        end

        # Anchor-drift weld: with a zero perturbation the stored walker
        # energy IS the raw accumulator, so a long incremental walk must
        # agree with a from-scratch recompute in the 1e-12 eV class
        weld_lat = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(6, 6, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:36]], adsorptions=:full)
        Random.seed!(97002)
        for i in 1:36
            weld_lat.components[1][i] = rand() < 0.5
        end
        weld_wk = LatticeWalker(weld_lat,
            energy=interacting_energy(weld_lat, inc_ham), iter=0)
        MC_grand_canonical_walk!(5000, weld_wk, inc_ham, 1.0e3, 0.0;
            p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=0.0,
            clusters_freq=2, swaps_freq=2, cluster_p=0.3, incremental=true)
        @test abs(ustrip(u"eV", weld_wk.energy -
                  interacting_energy(weld_wk.configuration, inc_ham))) <= 1e-12

        # Null-pair fixture at extreme coverage: a swap-only walk on a
        # nearly full lattice is dominated by equal-occupancy pairs, which
        # must contribute exactly zero delta (the weld exposes a corrupted
        # accumulator)
        null_lat = inc_lat()
        null_lat.components[1] .= true
        null_lat.components[1][3] = false
        Random.seed!(97003)
        null_wk = LatticeWalker(null_lat,
            energy=interacting_energy(null_lat, inc_ham), iter=0)
        MC_grand_canonical_walk!(2000, null_wk, inc_ham, 1.0e3, 0.0;
            p_move=1.0, p_insert=0.0, z0=1.0, energy_perturb=0.0,
            incremental=true)
        @test abs(ustrip(u"eV", null_wk.energy -
                  interacting_energy(null_wk.configuration, inc_ham))) <= 1e-12

        # Trait fallback: an unsupported Hamiltonian under incremental =
        # true recomputes fully each step and stays digit-identical to the
        # same-seed default
        cl_ham = ClusterLatticeHamiltonian(inc_ham,
            [ClusterInteraction(0.1u"eV", [(1, 2, 5)])])
        function inc_cl(inc)
            Random.seed!(97004)
            l = inc_lat()
            for i in 1:16
                l.components[1][i] = rand() < 0.5
            end
            w = LatticeWalker(l, energy=interacting_energy(l, cl_ham), iter=0)
            _, r, _, _, _, _ = MC_grand_canonical_walk!(500, w, cl_ham,
                1.0e3, 0.0; p_move=0.4, p_insert=0.3, z0=1.0,
                energy_perturb=1e-9, incremental=inc)
            return w.energy.val, sum(w.configuration.components[1]), r
        end
        eT, nT, rT = inc_cl(true)
        eF, nF, rF = inc_cl(false)
        @test eT == eF
        @test nT == nF
        @test rT == rF

        # Allocation ceiling on a compile-warmed incremental walk (generous
        # fixed byte bound, never time-based)
        alloc_lat = inc_lat()
        Random.seed!(97005)
        for i in 1:16
            alloc_lat.components[1][i] = rand() < 0.5
        end
        alloc_wk = LatticeWalker(alloc_lat,
            energy=interacting_energy(alloc_lat, inc_ham), iter=0)
        MC_grand_canonical_walk!(100, alloc_wk, inc_ham, 1.0e3, 0.0;
            p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
            incremental=true)
        alloc_bytes = @allocated MC_grand_canonical_walk!(100, alloc_wk,
            inc_ham, 1.0e3, 0.0; p_move=0.4, p_insert=0.3, z0=1.0,
            energy_perturb=1e-9, incremental=true)
        @test alloc_bytes < 2_000_000
    end

    @testset "driver controls (ideal-gas-referenced route)" begin
        dc_lat() = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1],
            components=[[false for _ in 1:16]], adsorptions=:full)
        dc_ham = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
        dc_save = SaveEveryN("t_dc_ig.csv", "t_dc_ig.traj", "t_dc_ig.ls",
                             1000000, 1000000, 1000000)
        dc_cleanup() = rm.(["t_dc_ig.csv", "t_dc_ig.traj", "t_dc_ig.ls"],
                           force=true)

        # Parameter-crafted stall: a near-zero reference fugacity
        # initializes every walker empty and a deletion-only move mix is
        # permanently guard-skipped
        function dc_stall(; kwargs...)
            Random.seed!(99003)
            ws = [LatticeWalker(deepcopy(dc_lat()), energy=0.0u"eV", iter=0)
                  for _ in 1:8]
            ls = LatticeGasWalkers(ws, dc_ham; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=10,
                reference_fugacity=1e-8, energy_perturbation=1e-9,
                allowed_fail_count=3)
            d, lsx, _ = ideal_gas_referenced_nested_sampling(ls, p,
                Int64(40), MCGrandCanonicalMoves(p_move=0.0, p_insert=0.0),
                dc_save; kwargs...)
            dc_cleanup()
            return d, lsx
        end
        d_stop, ls_stop = @test_logs (:warn, r"IG-ref GC-NS: Failed") match_mode=:any dc_stall(
            stop_on_stall=true)
        @test nrow(d_stop) == 0
        @test length(ls_stop.walkers) == 8
        d_def, _ = dc_stall()
        @test nrow(d_def) == 0

        # Tie-breaker scale warning: hand-computed bound on the 4x4
        # one-shell fixture is 0.04*16 + 0.01*64/2 = 0.96 in eV units; the
        # site-field wrapper adds its magnitude sum, the cluster wrapper
        # its coupling-times-embeddings sum, and an unknown Hamiltonian
        # type opts out (no warning)
        b0 = SamplingSchemes._perturbation_energy_bound(dc_ham, dc_lat())
        @test isapprox(b0, 0.96; rtol=1e-12)
        fld = collect(0.01 .* (1:16)) .* u"eV"
        bsf = SamplingSchemes._perturbation_energy_bound(
            SiteFieldLatticeHamiltonian(dc_ham, fld), dc_lat())
        @test isapprox(bsf, 0.96 + sum(0.01 .* (1:16)); rtol=1e-12)
        bcl = SamplingSchemes._perturbation_energy_bound(
            ClusterLatticeHamiltonian(dc_ham,
                [ClusterInteraction(0.1u"eV", [(1, 2, 5)])]), dc_lat())
        @test isapprox(bcl, 0.96 + 0.1; rtol=1e-12)
        @test SamplingSchemes._perturbation_energy_bound(:not_a_hamiltonian,
            dc_lat()) === nothing
        # Image-multiplicity cell (2x2, wrap-around bonds counted per image,
        # self-image entries included by the entries/2 form): each site
        # carries 4 ordered entries (two images per axis), 16 total, so
        # b = 0.04*4 + 0.01*16/2 = 0.24
        im_lat = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(2, 2, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1],
            components=[[false for _ in 1:4]], adsorptions=:full,
            image_multiplicity=true)
        @test isapprox(SamplingSchemes._perturbation_energy_bound(dc_ham,
            im_lat), 0.24; rtol=1e-12)
        # Multi-component Hamiltonian delegates to [1, 1], matching the
        # energy and delta paths
        @test isapprox(SamplingSchemes._perturbation_energy_bound(
            MLatticeHamiltonian(1, [dc_ham]), dc_lat()), 0.96; rtol=1e-12)

        # Below-floor perturbation warns at driver entry; an in-range one
        # never does
        function dc_short(delta)
            Random.seed!(99004)
            ws = [LatticeWalker(deepcopy(dc_lat()), energy=0.0u"eV", iter=0)
                  for _ in 1:10]
            ls = LatticeGasWalkers(ws, dc_ham; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=5,
                reference_fugacity=1.0, energy_perturbation=delta)
            d, _, _ = ideal_gas_referenced_nested_sampling(ls, p, Int64(3),
                MCGrandCanonicalMoves(), dc_save)
            dc_cleanup()
            return d
        end
        logs_lo, _ = Test.collect_test_logs() do
            dc_short(1e-15)
        end
        floor_warns = [l for l in logs_lo
                       if occursin("resolution floor", string(l.message))]
        @test length(floor_warns) == 1
        @test occursin("Recommended minimum", string(floor_warns[1].message))
        logs_hi, _ = Test.collect_test_logs() do
            dc_short(1e-9)
        end
        @test !any(occursin("resolution floor", string(l.message))
                   for l in logs_hi)
    end

    # ================================================================
    @testset "live-tail length guard" begin
        # Athermal fixture: 10-row ledger, K = 8, all energies zero, z = z0,
        # Skilling ω0 so Σω = 1 with the full tail. A truncated tail drops the
        # missing walkers' share of the residual mass X_n = (K/(K+1))^10.
        K_g, J_g = 8, 10
        df_g = DataFrame(iter=1:J_g, emax=zeros(J_g), num_particles=fill(3, J_g))
        w0_g = (K_g + 1) / K_g
        full_E, full_N = zeros(K_g), fill(3, K_g)
        half_E, half_N = zeros(4), fill(3, 4)

        # A K-length tail is warning-free
        s_full = @test_logs min_level = Base.CoreLogging.Warn gc_thermodynamic_stats_ideal_ref(
            df_g, 16, 1.0, [0.0], [300.0], K_g;
            ω0=w0_g, live_emax=full_E, live_numbers=full_N)
        # A truncated tail warns, naming both lengths and the tail-mass factor
        s_half = @test_logs (:warn, r"4 entries.*n_walkers = 8.*4/8") gc_thermodynamic_stats_ideal_ref(
            df_g, 16, 1.0, [0.0], [300.0], K_g;
            ω0=w0_g, live_emax=half_E, live_numbers=half_N)
        # The deficit the warning protects against, pinned to the closed form
        X_g = (K_g / (K_g + 1))^J_g
        @test isapprox(s_half.logXi[1, 1] - s_full.logXi[1, 1],
                       log(1 - X_g / 2); rtol=1e-12)
        # Dead-only calls are unaffected and warning-free
        s_dead = @test_logs min_level = Base.CoreLogging.Warn gc_thermodynamic_stats_ideal_ref(
            df_g, 16, 1.0, [0.0], [300.0], K_g)
        @test isfinite(s_dead.logXi[1, 1])
        # The effective-sample-size reduction's lattice method shares the guard
        @test_logs (:warn, r"4 entries.*n_walkers = 8") gc_effective_sample_size_ideal_ref(
            df_g, 16, 1.0, [0.0], [300.0], K_g;
            live_emax=half_E, live_numbers=half_N)
        e_full = @test_logs min_level = Base.CoreLogging.Warn gc_effective_sample_size_ideal_ref(
            df_g, 16, 1.0, [0.0], [300.0], K_g;
            live_emax=full_E, live_numbers=full_N)
        @test all(isfinite, e_full)
    end
end
