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

    function igref_run(template, ham, z0, n_walkers, n_steps; seed=7, mc_steps=100)
        Random.seed!(seed)
        walkers = [LatticeWalker(deepcopy(template), energy=0.0u"eV", iter=0)
                   for _ in 1:n_walkers]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=mc_steps, reference_fugacity=z0)
        mc_routine = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)
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

        # z0 must be positive
        walker = LatticeWalker(deepcopy(igref_template), energy=0.0u"eV", iter=0)
        @test_throws ArgumentError MC_grand_canonical_walk!(
            1, walker, ham0, 1e3, 0.0; z0=0.0)
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
        stats_live = gc_thermodynamic_stats_ideal_ref(
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
end
