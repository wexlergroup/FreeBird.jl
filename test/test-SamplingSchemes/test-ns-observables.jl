@testset "NS observables" begin
    using Random

    # Shared single-component square-lattice builder
    function obs_square_lattice(d1, d2)
        MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(d1, d2, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:d1*d2]],
            adsorptions=:full)
    end

    # ================================================================
    @testset "order_parameter_c2x2" begin
        lat = obs_square_lattice(4, 4)

        # Empty and full lattices carry no c(2x2) order
        @test order_parameter_c2x2(lat) == 0.0
        lat.components[1] .= true
        @test order_parameter_c2x2(lat) == 0.0

        # Single particle: |±1| / M
        lat.components[1] .= false
        lat.components[1][1] = true
        @test order_parameter_c2x2(lat) == 1 / 16

        # Perfect checkerboards (both sublattices) at half filling: 1/2.
        # Site ordering: dimension 1 fastest, so i0 = (s-1) % d1, j0 = (s-1) ÷ d1.
        checker = [iseven(((s - 1) % 4) + ((s - 1) ÷ 4)) for s in 1:16]
        lat.components[1] .= checker
        @test order_parameter_c2x2(lat) == 0.5
        lat.components[1] .= .!checker
        @test order_parameter_c2x2(lat) == 0.5

        # A column stripe alternates parity along the column and cancels
        lat.components[1] .= [(s - 1) % 4 == 0 for s in 1:16]
        @test order_parameter_c2x2(lat) == 0.0

        # Odd in-plane dimensions: the sublattices do not tile evenly
        @test_throws ArgumentError order_parameter_c2x2(obs_square_lattice(3, 4))
        @test_throws ArgumentError order_parameter_c2x2(obs_square_lattice(4, 3))

        # Three-dimensional supercell
        lat3d = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 2),
            periodicity=(true, true, true),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:32]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_c2x2(lat3d)

        # Multi-site basis
        lat2b = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)],
            supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[0.8],
            components=[[false for _ in 1:32]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_c2x2(lat2b)
    end

    # ================================================================
    obs_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
    obs_save = SaveEveryN("t_obs.csv", "t_obs.traj", "t_obs.ls",
                          1000000, 1000000, 1000000)
    obs_cleanup() = rm.(["t_obs.csv", "t_obs.traj", "t_obs.ls"], force=true)
    count_N(cfg) = Float64(sum(cfg.components[1]))

    @testset "observable hook: validation" begin
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:8]
        ls = LatticeGasWalkers(walkers, obs_ham)
        params = LatticeNestedSamplingParameters(mc_steps=10, energy_perturbation=1e-9)

        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=[:a => order_parameter_c2x2, :a => order_parameter_c2x2])
        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=[:emax => order_parameter_c2x2])
        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=Pair{Symbol,Function}[])
        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=[:bad => (cfg -> "not a number")])
        # Non-finite probe values warn but do not throw
        @test_logs (:warn, r"non-finite") min_level=Base.CoreLogging.Warn match_mode=:any nested_sampling(
            ls, params, Int64(2), MCRandomWalkClone(), obs_save;
            observables=[:nanobs => (cfg -> NaN)])
        obs_cleanup()
    end

    # ================================================================
    @testset "observable hook: exact dead-point pairing (igref route)" begin
        Random.seed!(42)
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:30]
        ls = LatticeGasWalkers(walkers, obs_ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=30, reference_fugacity=1.0)
        mc = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)

        df, final_ls, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(400), mc, obs_save;
            observables=[:n_check => count_N, :psi => order_parameter_c2x2])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "num_particles", "n_check", "psi"]
        @test eltype(df.n_check) == Float64
        @test eltype(df.psi) == Float64
        # The load-bearing check: the observable column, evaluated on the
        # captured culled walker, reproduces the independently recorded
        # num_particles column row by row over a full run
        @test df.n_check == Float64.(df.num_particles)
        @test all(0.0 .<= df.psi .<= 0.5)

        # Schema is unchanged when observables are not requested
        Random.seed!(42)
        walkers2 = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:30]
        ls2 = LatticeGasWalkers(walkers2, obs_ham; assign_energy=false)
        params2 = IdealGasReferencedGCNSParameters(mc_steps=30, reference_fugacity=1.0)
        df2, _, _ = ideal_gas_referenced_nested_sampling(
            ls2, params2, Int64(50), mc, obs_save)
        obs_cleanup()
        @test names(df2) == ["iter", "emax", "num_particles"]
    end

    # ================================================================
    @testset "observable hook: exact dead-point pairing (omega-sorted route)" begin
        Random.seed!(11)
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:20]
        ls = LatticeGasWalkers(walkers, obs_ham; assign_energy=false)
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=30, chemical_potential=-0.02)
        mc = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)

        df, _, _ = grand_canonical_nested_sampling(
            ls, gc_params, Int64(200), mc, obs_save;
            observables=[:n_check => count_N])
        obs_cleanup()

        @test names(df) == ["iter", "omega", "energy", "num_particles", "n_check"]
        @test df.n_check == Float64.(df.num_particles)
    end

    # ================================================================
    @testset "observable hook: exact dead-point pairing (canonical route)" begin
        Random.seed!(7)
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:20]
        # Fixed-N ladder at N = 5: place 5 particles per walker
        for w in walkers
            occ = vcat(fill(true, 5), fill(false, 11))
            shuffle!(occ)
            w.configuration.components[1] .= occ
        end
        ls = LatticeGasWalkers(walkers, obs_ham; perturb_energy=1e-9)
        params = LatticeNestedSamplingParameters(mc_steps=30,
            energy_perturbation=1e-9, allowed_fail_count=100000)

        df, _, _ = nested_sampling(ls, params, Int64(200),
            MCRandomWalkClone(), obs_save;
            observables=[:e_check => (cfg -> interacting_energy(cfg, obs_ham).val),
                         :psi => order_parameter_c2x2])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "e_check", "psi"]
        # N is conserved on the canonical route, so pair the recomputed
        # (unperturbed) energy against the recorded emax instead; they differ
        # only by the tie-breaking perturbation
        @test all(isapprox.(df.e_check, df.emax; atol=1e-8))
        @test all(0.0 .<= df.psi .<= 0.5)
    end

    # ================================================================
    @testset "igref stats: var_U, cov_UN, and observable averages" begin
        kb = 8.617333262e-5
        df = DataFrame(iter=[1, 2, 3], emax=[0.5, 0.3, 0.1],
                       num_particles=[2, 1, 0], a=[1.0, 2.0, 4.0])
        live_E = [0.05, 0.02]
        live_N = [1, 2]
        live_a = [8.0, 16.0]
        K, z0, M = 4, 1.0, 16
        ω0 = (K + 1) / K
        μs = [-0.02, 0.03]
        Ts = [300.0, 600.0]

        stats = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K;
            ω0=ω0, live_emax=live_E, live_numbers=live_N,
            observable_cols=[:a], live_observables=Dict(:a => live_a))

        # Independent linear-space reference for the closed-form weighted sums
        w0 = ω0 * (1 / (K + 1)) .* (K / (K + 1)) .^ df.iter
        wt = fill((K / (K + 1))^3 / K, 2)      # tail carries no ω0 factor
        w_all = vcat(w0, wt)
        E_all = vcat(df.emax, live_E)
        N_all = vcat(Float64.(df.num_particles), Float64.(live_N))
        a_all = vcat(df.a, live_a)
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            β = 1 / (kb * T)
            w = w_all .* exp.(β * μ .* N_all .- β .* E_all)   # z0 = 1
            sw = sum(w)
            u = sum(w .* E_all) / sw
            n = sum(w .* N_all) / sw
            @test stats.mean_U[i, j] ≈ u rtol = 1e-12
            @test stats.var_U[i, j] ≈ sum(w .* E_all .^ 2) / sw - u^2 rtol = 1e-10
            @test stats.cov_UN[i, j] ≈ sum(w .* E_all .* N_all) / sw - u * n rtol = 1e-10
            @test stats.observables[:a][i, j] ≈ sum(w .* a_all) / sw rtol = 1e-12
        end

        # Self-consistency: averaging the ledger's own columns reproduces
        # mean_N and mean_U
        stats2 = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K;
            ω0=ω0, live_emax=live_E, live_numbers=live_N,
            observable_cols=[:num_particles, :emax],
            live_observables=Dict(:num_particles => Float64.(live_N),
                                  :emax => live_E))
        @test stats2.observables[:num_particles] ≈ stats2.mean_N rtol = 1e-12
        @test stats2.observables[:emax] ≈ stats2.mean_U rtol = 1e-12

        # Backward compatibility: pre-existing fields keep their order, new
        # fields are appended, positional destructuring still works
        stats3 = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K)
        @test isempty(stats3.observables)
        @test keys(stats3) == (:logXi, :mean_N, :var_N, :mean_U, :N_eff,
                               :var_U, :cov_UN, :observables)
        lX, mN, vN, mU, Ne = stats3
        @test lX == stats3.logXi && Ne == stats3.N_eff

        # Validation traps
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:nope])
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a, :a],
            live_emax=live_E, live_numbers=live_N,
            live_observables=Dict(:a => live_a))
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a],
            live_emax=live_E, live_numbers=live_N)
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a],
            live_emax=live_E, live_numbers=live_N,
            live_observables=Dict(:b => live_a))
        @test_throws DimensionMismatch gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a],
            live_emax=live_E, live_numbers=live_N,
            live_observables=Dict(:a => [1.0]))
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K;
            live_observables=Dict(:a => live_a))
    end

    # ================================================================
    @testset "fixed-N stats: var_U, cov_UN, and observable averages" begin
        kb = 8.617333262e-5
        M = 4
        K = 3
        ω0 = (K + 1) / K
        N_values = [0, 1, 2]
        T_grid = [300.0, 600.0] .* u"K"
        μ_grid = [-0.05, 0.02] .* u"eV"

        df0 = DataFrame(iter=Int[], emax=Float64[])
        df1 = DataFrame(iter=[1, 2], emax=[0.3, 0.1], psi=[0.2, 0.4])
        df2 = DataFrame(iter=Int[], emax=Float64[])   # single-config convention
        live_E = [Float64[], [0.05, 0.01], [0.02, 0.02, 0.02]]
        live_psi = [Dict(:psi => [0.0]),
                    Dict(:psi => [0.6, 0.8]),
                    Dict(:psi => [0.5, 0.5, 0.5])]

        stats = gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, ω0=ω0, live_emax=live_E,
            observable_cols=[:psi], live_observables=live_psi)

        # Independent linear-space reference. Sector 1: two dead points plus
        # a two-entry tail (n_iters = 2). Sector 2: empty ladder, tail mass
        # exactly 1 split over three copies of the sector energy.
        w1 = vcat(ω0 * (1 / (K + 1)) .* (K / (K + 1)) .^ df1.iter,
                  fill(ω0 * (K / (K + 1))^2 / 2, 2))
        E1 = vcat(df1.emax, live_E[2])
        a1 = vcat(df1.psi, live_psi[2][:psi])
        binom = [1.0, 4.0, 6.0]                       # C(4, N)
        for (j, T) in enumerate([300.0, 600.0]), (k, μ) in enumerate([-0.05, 0.02])
            β = 1 / (kb * T)
            b1 = w1 .* exp.(-β .* E1)
            z1 = sum(b1)
            e1 = sum(b1 .* E1) / z1
            e1sq = sum(b1 .* E1 .^ 2) / z1
            p1 = sum(b1 .* a1) / z1
            z2 = exp(-β * 0.02)
            e2, e2sq, p2 = 0.02, 0.02^2, 0.5
            wN = binom .* [1.0, z1, z2] .* exp.(β * μ .* [0.0, 1.0, 2.0])
            sw = sum(wN)
            U = (wN[2] * e1 + wN[3] * e2) / sw
            Nbar = (wN[2] * 1 + wN[3] * 2) / sw
            @test stats.logXi[k, j] ≈ log(sw) rtol = 1e-10
            @test stats.mean_N[k, j] ≈ Nbar rtol = 1e-10
            @test stats.mean_U[k, j] ≈ U rtol = 1e-10
            @test stats.var_U[k, j] ≈
                  (wN[2] * e1sq + wN[3] * e2sq) / sw - U^2 rtol = 1e-8
            @test stats.cov_UN[k, j] ≈
                  (wN[2] * 1 * e1 + wN[3] * 2 * e2) / sw - U * Nbar rtol = 1e-8
            # The N = 0 sector contributes psi = 0 with weight wN[1]
            @test stats.observables[:psi][k, j] ≈
                  (wN[2] * p1 + wN[3] * p2) / sw rtol = 1e-10
        end

        # Backward compatibility: field order preserved, new fields appended
        stats_nc = gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, ω0=ω0, live_emax=live_E)
        @test keys(stats_nc) == (:logXi, :mean_N, :var_N, :mean_U,
                                 :log_Z_N, :N_values, :var_U, :cov_UN,
                                 :observables)
        @test isempty(stats_nc.observables)
        @test stats_nc.logXi ≈ stats.logXi rtol = 1e-14

        # Validation traps
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi])
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, observable_cols=[:psi], live_observables=live_psi)
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=live_psi[1:2])
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=[Dict(:other => [0.0]), live_psi[2], live_psi[3]])
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=[live_psi[1], Dict(:psi => [0.6]), live_psi[3]])
        df1_nocol = DataFrame(iter=[1, 2], emax=[0.3, 0.1])
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1_nocol, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=live_psi)
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, live_observables=live_psi)
    end
end
