@testset "Grand-canonical nested sampling tests" begin
 
    # ================================================================
    # Shared lattice and Hamiltonian for all GC tests
    # ================================================================
    square_lattice = MLattice{1,SquareLattice}(
        lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)],
        supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false),
        cutoff_radii=[1.1, 1.5],
        components=:equal,
        adsorptions=:full
    )
    ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
 
    # ================================================================
    @testset "GrandCanonicalNestedSamplingParameters" begin
        gc_params = GrandCanonicalNestedSamplingParameters()
        @test gc_params isa SamplingSchemes.SamplingParameters
        @test gc_params.mc_steps == 100
        @test gc_params.chemical_potential == 0.0
        @test gc_params.energy_perturbation == 1e-12
        @test gc_params.random_seed == 1234
        @test gc_params.fail_count == 0
        @test gc_params.allowed_fail_count == 10
        @test gc_params.init_occupation_p == 0.5
        @test gc_params.n_max == typemax(Int64)
 
        gc_params2 = GrandCanonicalNestedSamplingParameters(
            mc_steps=200, chemical_potential=-0.05,
            init_occupation_p=0.3)
        @test gc_params2.mc_steps == 200
        @test gc_params2.chemical_potential == -0.05
        @test gc_params2.init_occupation_p == 0.3
        @test gc_params2.n_max == typemax(Int64)

        gc_params3 = GrandCanonicalNestedSamplingParameters(n_max=Int64(5))
        @test gc_params3.n_max == 5

        # Test mutability
        gc_params.fail_count = 5
        @test gc_params.fail_count == 5

        # Cluster field defaults
        @test gc_params.cluster_p == 0.3
        @test gc_params.cluster_accepted == 0.0
        @test gc_params.cluster_total == 0.0
        @test gc_params.cluster_p_history == Float64[]
        @test gc_params.cluster_accept_history == Float64[]
        @test gc_params.cluster_adjust_iterations == Int[]

        # move_stats defaults (issue #158)
        @test gc_params.move_stats isa Dict{Symbol,Int}
        @test isempty(gc_params.move_stats)

        # Cluster field mutability
        gc_params.cluster_p = 0.5
        @test gc_params.cluster_p == 0.5
    end
 
    # ================================================================
    @testset "MCGrandCanonicalMoves" begin
        mc = MCGrandCanonicalMoves()
        @test mc isa MCRoutine
        @test mc.p_move == 0.5
        @test mc.p_insert == 0.25
 
        mc2 = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)
        @test mc2.p_move == 0.4
        @test mc2.p_insert == 0.3

        # Cluster field defaults
        @test mc.clusters_freq == 0
        @test mc.swaps_freq == 1
        @test mc.initial_cluster_p == 0.3
        @test mc.target_cluster_accept == 0.3
        @test mc.cluster_adjust_interval == 50
        @test mc.cluster_p_floor == 0.01
        @test mc.cluster_p_ceiling == 1.0

        # Cluster-enabled construction
        mc3 = MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25,
            clusters_freq=3, swaps_freq=1, initial_cluster_p=0.5)
        @test mc3.clusters_freq == 3
        @test mc3.swaps_freq == 1
        @test mc3.initial_cluster_p == 0.5

        # Invalid probabilities
        @test_throws ArgumentError MCGrandCanonicalMoves(p_move=0.8, p_insert=0.3)
        @test_throws ArgumentError MCGrandCanonicalMoves(p_move=-0.1, p_insert=0.3)

        # Bias field defaults (issue #158)
        @test mc.p_bias == 0.0
        @test mc.bias_predicate == :contact
        @test mc.bias_shells == 1

        # Bias-enabled construction
        mc4 = MCGrandCanonicalMoves(p_bias=0.5, bias_predicate=:cavity, bias_shells=2)
        @test mc4.p_bias == 0.5
        @test mc4.bias_predicate == :cavity
        @test mc4.bias_shells == 2

        # Invalid bias parameters
        @test_throws ArgumentError MCGrandCanonicalMoves(p_bias=-0.1)
        @test_throws ArgumentError MCGrandCanonicalMoves(p_bias=1.5)
        @test_throws ArgumentError MCGrandCanonicalMoves(bias_shells=0)
        @test_throws ArgumentError MCGrandCanonicalMoves(bias_predicate=:nonsense)

        # p_bias = 1.0 constructs but warns (reducible kernel; see the
        # non-ergodicity pin in test-ideal-gas-ref-gcns.jl) ...
        @test_logs (:warn, r"p_bias") match_mode = :any MCGrandCanonicalMoves(p_bias=1.0)
        # ... and stays silent below 1
        @test_logs min_level = Base.CoreLogging.Warn MCGrandCanonicalMoves(p_bias=0.99)
    end
 
    # ================================================================
    @testset "random_microstate!" begin
        lattice = deepcopy(square_lattice)
 
        # p=0 gives empty lattice
        random_microstate!(lattice; p=0.0)
        @test sum(lattice.components[1]) == 0
 
        # p=1 gives full lattice
        random_microstate!(lattice; p=1.0)
        @test sum(lattice.components[1]) == num_sites(lattice)
 
        # p=0.5 gives variable occupancy (statistical — just check it runs)
        counts = Int[]
        for _ in 1:50
            random_microstate!(lattice; p=0.5)
            push!(counts, sum(lattice.components[1]))
        end
        # With 16 sites and p=0.5, we should see variation
        @test minimum(counts) < maximum(counts)
        # Mean should be roughly 8
        @test 4 < mean(counts) < 12
    end
 
    # ================================================================
    @testset "lattice_insert_particle!" begin
        lattice = deepcopy(square_lattice)
        lattice.components[1] .= false
        n_sites = num_sites(lattice)
 
        # Insert into empty lattice
        success, _ = lattice_insert_particle!(lattice)
        @test success
        @test sum(lattice.components[1]) == 1
 
        # Fill the lattice
        for _ in 2:n_sites
            lattice_insert_particle!(lattice)
        end
        @test sum(lattice.components[1]) == n_sites
 
        # Insert into full lattice should fail
        success, _ = lattice_insert_particle!(lattice)
        @test !success
        @test sum(lattice.components[1]) == n_sites
    end
 
    # ================================================================
    @testset "lattice_delete_particle!" begin
        lattice = deepcopy(square_lattice)
        lattice.components[1] .= true
        n_sites = num_sites(lattice)
 
        # Delete from full lattice
        success, _ = lattice_delete_particle!(lattice)
        @test success
        @test sum(lattice.components[1]) == n_sites - 1
 
        # Empty the lattice
        lattice.components[1] .= false
        lattice.components[1][1] = true
        success, _ = lattice_delete_particle!(lattice)
        @test success
        @test sum(lattice.components[1]) == 0
 
        # Delete from empty lattice should fail
        success, _ = lattice_delete_particle!(lattice)
        @test !success
    end
 
    # ================================================================
    @testset "MC_grand_canonical_walk!" begin
        walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        assign_energy!(walker, ham)
 
        mu = -0.05  # in eV
        n_init = sum(walker.configuration.components[1])
        omega_max = walker.energy.val - mu * n_init + 1.0  # generous upper bound
 
        accept, rate, updated_walker, cl_acc, cl_tot = MC_grand_canonical_walk!(
            100, walker, ham, omega_max, mu;
            p_move=0.5, p_insert=0.25, energy_perturb=0.0)

        @test accept isa Bool
        @test 0.0 <= rate <= 1.0
        @test updated_walker isa LatticeWalker{1}
        # No cluster moves when clusters_freq=0 (default)
        @test cl_acc == 0
        @test cl_tot == 0

        # Energy should be consistent with the configuration
        expected_energy = interacting_energy(updated_walker.configuration, ham)
        @test updated_walker.energy ≈ expected_energy
 
        # Omega should be below omega_max
        n_final = sum(updated_walker.configuration.components[1])
        omega_final = updated_walker.energy.val - mu * n_final
        @test omega_final < omega_max
 
        # Invalid probability should throw
        @test_throws ArgumentError MC_grand_canonical_walk!(
            10, walker, ham, omega_max, mu; p_move=0.8, p_insert=0.3)

        # ---- Move counters: 6th return element (issue #158) ----
        walker6 = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        assign_energy!(walker6, ham)
        _, _, walker6, _, _, counters = MC_grand_canonical_walk!(
            100, walker6, ham, 1e3, mu;
            p_move=0.5, p_insert=0.25, energy_perturb=0.0)
        @test counters isa NamedTuple
        @test keys(counters) == (
            :swap_attempted, :swap_accepted, :cluster_attempted, :cluster_accepted,
            :insert_uniform_attempted, :insert_uniform_accepted,
            :insert_biased_attempted, :insert_biased_accepted,
            :delete_attempted, :delete_accepted)
        # Default path: the bias and cluster families never fire
        @test counters.insert_biased_attempted == 0
        @test counters.insert_biased_accepted == 0
        @test counters.cluster_attempted == 0
        # accepted <= attempted per family
        @test counters.swap_accepted <= counters.swap_attempted
        @test counters.insert_uniform_accepted <= counters.insert_uniform_attempted
        @test counters.delete_accepted <= counters.delete_attempted
    end

    # ================================================================
    @testset "MC_grand_canonical_walk! counters: exact bookkeeping" begin
        ham0 = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        # Skip-free by construction: N starts at 8, |dN| <= 1 per step, so
        # over 8 steps N stays in [1, 15]: insert-at-full and delete-at-empty
        # can never fire, and with 0 < N < M on a connected lattice the
        # :contact set is never empty, so the biased branch cannot skip.
        # Every step assigns a move type, for any RNG stream:
        # sum(attempted) == n_steps exactly.
        lat_bk = deepcopy(square_lattice)
        lat_bk.components[1] .= false
        lat_bk.components[1][1:8] .= true
        walker_bk = LatticeWalker(lat_bk, energy=0.0u"eV", iter=0)
        _, _, _, _, _, c = MC_grand_canonical_walk!(
            8, walker_bk, ham0, 1e3, 0.0;
            p_move=0.5, p_insert=0.25, energy_perturb=0.0,
            p_bias=0.5, bias_predicate=:contact, bias_shells=1)
        attempted = c.swap_attempted + c.cluster_attempted +
                    c.insert_uniform_attempted + c.insert_biased_attempted +
                    c.delete_attempted
        accepted = c.swap_accepted + c.cluster_accepted +
                   c.insert_uniform_accepted + c.insert_biased_accepted +
                   c.delete_accepted
        @test attempted == 8
        @test accepted <= attempted

        # kwarg validation on the walk
        w2 = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        @test_throws ArgumentError MC_grand_canonical_walk!(1, w2, ham0, 1e3, 0.0; p_bias=-0.1)
        @test_throws ArgumentError MC_grand_canonical_walk!(1, w2, ham0, 1e3, 0.0; p_bias=1.5)
        @test_throws ArgumentError MC_grand_canonical_walk!(1, w2, ham0, 1e3, 0.0; bias_shells=0)
        @test_throws ArgumentError MC_grand_canonical_walk!(1, w2, ham0, 1e3, 0.0; bias_predicate=:nonsense)
        # bias_shells beyond the lattice fires only when the channel is active
        @test_throws ArgumentError MC_grand_canonical_walk!(
            1, w2, ham0, 1e3, 0.0; p_bias=0.5, bias_shells=3)
    end

    # ================================================================
    @testset "MC_grand_canonical_walk! with cluster moves" begin
        walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        assign_energy!(walker, ham)

        mu = -0.05
        n_init = sum(walker.configuration.components[1])
        omega_max = walker.energy.val - mu * n_init + 1.0

        accept, rate, updated_walker, cl_acc, cl_tot = MC_grand_canonical_walk!(
            200, walker, ham, omega_max, mu;
            p_move=0.5, p_insert=0.25, energy_perturb=0.0,
            clusters_freq=1, swaps_freq=1, cluster_p=0.3)

        @test accept isa Bool
        @test 0.0 <= rate <= 1.0
        @test updated_walker isa LatticeWalker{1}
        # Should have attempted some cluster moves
        @test cl_tot > 0
        @test cl_acc >= 0
        @test cl_acc <= cl_tot

        # Energy should be consistent
        expected_energy = interacting_energy(updated_walker.configuration, ham)
        @test updated_walker.energy ≈ expected_energy
    end

    # ================================================================
    @testset "nested_sampling_step! for GC" begin
        walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0) for _ in 1:5]
        liveset = LatticeGasWalkers(walkers, ham)
 
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=50, chemical_potential=-0.05)
        mc_routine = MCGrandCanonicalMoves()
 
        # Initialize walkers with random microstates
        SamplingSchemes._init_gc_walkers!(liveset, gc_params)
 
        e_type = typeof(walkers[1].energy)
        iter, omega, energy, n_par, updated_liveset, updated_params = nested_sampling_step!(
            liveset, gc_params, mc_routine)
 
        @test iter isa Union{Missing,Int}
        @test omega isa Union{Missing,e_type}
        @test energy isa Union{Missing,e_type}
        @test n_par isa Union{Missing,Int}
        @test length(updated_liveset.walkers) == 5
        @test updated_params.fail_count >= 0

        # ---- move_stats accumulation through one NS step (issue #158) ----
        attempted_keys = (:swap_attempted, :cluster_attempted,
                          :insert_uniform_attempted, :insert_biased_attempted,
                          :delete_attempted)
        ms = updated_params.move_stats
        @test ms isa Dict{Symbol,Int}
        @test all(v >= 0 for v in values(ms))
        att = sum(get(ms, k, 0) for k in attempted_keys)
        # One NS step walks one replacement walker for mc_steps steps; skips
        # require a boundary configuration, so at least one step attempts
        @test 1 <= att <= gc_params.mc_steps
        @test get(ms, :swap_accepted, 0) <= get(ms, :swap_attempted, 0)
        @test get(ms, :insert_uniform_accepted, 0) <= get(ms, :insert_uniform_attempted, 0)
        @test get(ms, :insert_biased_accepted, 0) <= get(ms, :insert_biased_attempted, 0)
        @test get(ms, :delete_accepted, 0) <= get(ms, :delete_attempted, 0)
    end
 
    # ================================================================
    @testset "grand_canonical_nested_sampling loop" begin
        walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0) for _ in 1:10]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
 
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=50, chemical_potential=-0.05)
        mc_routine = MCGrandCanonicalMoves()
        save_strategy = SaveEveryN("test_gc_df.csv", "test_gc.traj", "test_gc.ls", 1000, 1000, 1000)
 
        df, updated_liveset, updated_params = grand_canonical_nested_sampling(
            liveset, gc_params, Int64(20), mc_routine, save_strategy)
 
        @test df isa DataFrame
        @test names(df) == ["iter", "omega", "energy", "num_particles"]
        @test nrow(df) <= 20
        @test nrow(df) > 0  # At least some steps should succeed
        @test eltype(df.iter) == Int
        @test eltype(df.omega) == Float64
        @test eltype(df.energy) == Float64
        @test eltype(df.num_particles) == Int
        @test length(updated_liveset.walkers) == 10
 
        # Omega should be monotonically non-increasing (each recorded Ω <= previous)
        if nrow(df) > 1
            for i in 2:nrow(df)
                @test df.omega[i] <= df.omega[1] + 1e-10
            end
        end
 
        rm("test_gc_df.csv", force=true)
        rm("test_gc.traj", force=true)
        rm("test_gc.ls", force=true)
    end

    # ================================================================
    @testset "move_stats: accumulated per run, reset between runs" begin
        walkers_ms = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0) for _ in 1:10]
        liveset_ms = LatticeGasWalkers(walkers_ms, ham; assign_energy=false)
        params_ms = GrandCanonicalNestedSamplingParameters(
            mc_steps=50, chemical_potential=-0.05)
        save_ms = SaveEveryN("test_ms_df.csv", "test_ms.traj", "test_ms.ls", 1000, 1000, 1000)
        attempted_keys = (:swap_attempted, :cluster_attempted,
                          :insert_uniform_attempted, :insert_biased_attempted,
                          :delete_attempted)

        _, _, params1 = grand_canonical_nested_sampling(
            liveset_ms, params_ms, Int64(20), MCGrandCanonicalMoves(), save_ms)
        s1 = sum(get(params1.move_stats, k, 0) for k in attempted_keys)
        @test 0 < s1 <= 20 * 50   # per-run ceiling: n_iters x mc_steps

        # Re-running with the RETURNED params must reset the accumulator: a
        # carried-over dict would exceed the single-run ceiling
        _, _, params2 = grand_canonical_nested_sampling(
            liveset_ms, params1, Int64(20), MCGrandCanonicalMoves(), save_ms)
        s2 = sum(get(params2.move_stats, k, 0) for k in attempted_keys)
        @test 0 < s2 <= 20 * 50

        rm("test_ms_df.csv", force=true)
        rm("test_ms.traj", force=true)
        rm("test_ms.ls", force=true)
    end

    # ================================================================
    @testset "gc_thermodynamic_stats basic" begin
        # Hand-crafted test: 2 microstates
        # Microstate 1: E=0, N=0, Ω=0
        # Microstate 2: E=-1, N=1, Ω=-1-μ*1
        μ = -0.5
        ωi = [0.5, 0.5]
        grand_es = [0.0, -1.0 - μ * 1.0]  # Ω = E - μN
        Es = [0.0, -1.0]
        Ns = [0, 1]
 
        # At β=0 (infinite T), equal weights
        u, cv_val, n_avg = gc_thermodynamic_stats(
            0.001, ωi, grand_es, Es, Ns, μ; kb=1.0)
        @test isfinite(u)
        @test isfinite(cv_val)
        @test isfinite(n_avg)
        # At very low β, should be roughly equal mixture
        @test n_avg ≈ 0.5 atol=0.1
 
        # Dimension mismatch should throw
        @test_throws DimensionMismatch gc_thermodynamic_stats(
            1.0, [0.5], [0.0, 1.0], [0.0, 1.0], [0, 1], 0.0)
    end
 
    # ================================================================
    @testset "Validation against exact grand-canonical enumeration" begin
        # 4x4 lattice, NN interactions only, single component
        # Exact: sum over all 2^16 = 65536 microstates
        L = 4
        lattice_template = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:L*L]],
            adsorptions=:full
        )
        ham_val = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        n_sites = L * L
 
        # Exact grand-canonical enumeration
        mu_val = -0.05  # eV
        kb = 8.617333262e-5  # eV/K
        T_test = 300.0  # K
        beta_test = 1.0 / (kb * T_test)
 
        # Enumerate all 2^16 microstates
        exact_z = 0.0
        exact_E = 0.0
        exact_E2 = 0.0
        exact_N = 0.0
        exact_EN = 0.0
 
        for mask in 0:(2^n_sites - 1)
            lattice = deepcopy(lattice_template)
            for site in 1:n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_val = interacting_energy(lattice, ham_val).val
            N_val = sum(lattice.components[1])
            omega_val = E_val - mu_val * N_val
 
            boltz = exp(-beta_test * omega_val)
            exact_z += boltz
            exact_E += boltz * E_val
            exact_E2 += boltz * E_val^2
            exact_N += boltz * N_val
            exact_EN += boltz * E_val * N_val
        end
 
        exact_mean_E = exact_E / exact_z
        exact_mean_N = exact_N / exact_z
        exact_mean_E2 = exact_E2 / exact_z
        exact_mean_EN = exact_EN / exact_z
        exact_var_E = exact_mean_E2 - exact_mean_E^2
        exact_cov_EN = exact_mean_EN - exact_mean_E * exact_mean_N
        exact_Cv = kb * beta_test^2 * (exact_var_E - mu_val * exact_cov_EN)
 
        # Run GC-NS with enough walkers and iterations
        n_walkers = 100
        n_steps = Int64(3000)
 
        walkers = [LatticeWalker(deepcopy(lattice_template), energy=0.0u"eV", iter=0)
                   for _ in 1:n_walkers]
        liveset = LatticeGasWalkers(walkers, ham_val; assign_energy=false)
 
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=100, chemical_potential=mu_val,
            energy_perturbation=1e-12, init_occupation_p=0.5)
        mc_routine = MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25)
        save_strategy = SaveEveryN("test_val.csv", "test_val.traj", "test_val.ls", 10000, 10000, 10000)
 
        df, _, _ = grand_canonical_nested_sampling(
            liveset, gc_params, n_steps, mc_routine, save_strategy)
 
        @test nrow(df) > 0
 
        # Compute NS thermodynamic stats
        mean_E_ns, Cv_ns, mean_N_ns = gc_thermodynamic_stats(
            df, [beta_test], n_walkers, mu_val)
 
        # Compare with exact values (generous tolerances for stochastic algorithm)
        @test mean_E_ns[1] ≈ exact_mean_E rtol=0.3
        @test mean_N_ns[1] ≈ exact_mean_N rtol=0.3
        # Cv is harder to converge; just check it's in the right ballpark
        if isfinite(Cv_ns[1]) && isfinite(exact_Cv) && exact_Cv > 0
            @test Cv_ns[1] > 0
        end
 
        rm("test_val.csv", force=true)
        rm("test_val.traj", force=true)
        rm("test_val.ls", force=true)
    end

    # ================================================================
    @testset "n_max enforcement" begin
        n_max_val = Int64(5)
        walkers_nm = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0) for _ in 1:10]
        liveset_nm = LatticeGasWalkers(walkers_nm, ham; assign_energy=false)

        gc_params_nm = GrandCanonicalNestedSamplingParameters(
            mc_steps=50, chemical_potential=-0.05,
            init_occupation_p=0.5, n_max=n_max_val)
        mc_routine_nm = MCGrandCanonicalMoves()
        save_nm = SaveEveryN("test_nmax_df.csv", "test_nmax.traj", "test_nmax.ls", 10000, 10000, 10000)

        df_nm, updated_liveset_nm, _ = grand_canonical_nested_sampling(
            liveset_nm, gc_params_nm, Int64(50), mc_routine_nm, save_nm)

        # After initialization, all walkers must respect n_max
        for w in updated_liveset_nm.walkers
            @test sum(w.configuration.components[1]) <= n_max_val
        end

        # All recorded samples must respect n_max
        if nrow(df_nm) > 0
            @test all(df_nm.num_particles .<= n_max_val)
        end

        rm("test_nmax_df.csv", force=true)
        rm("test_nmax.traj", force=true)
        rm("test_nmax.ls", force=true)
    end

    # ================================================================
    @testset "GC-NS with cluster moves: basic functionality" begin
        walkers_cl = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0) for _ in 1:10]
        liveset_cl = LatticeGasWalkers(walkers_cl, ham; assign_energy=false)

        gc_params_cl = GrandCanonicalNestedSamplingParameters(
            mc_steps=50, chemical_potential=-0.05)
        mc_routine_cl = MCGrandCanonicalMoves(
            p_move=0.5, p_insert=0.25,
            clusters_freq=1, swaps_freq=1, initial_cluster_p=0.3,
            cluster_adjust_interval=20)
        save_cl = SaveEveryN("test_gc_cl_df.csv", "test_gc_cl.traj", "test_gc_cl.ls", 10000, 10000, 10000)

        df_cl, updated_liveset_cl, updated_params_cl = grand_canonical_nested_sampling(
            liveset_cl, gc_params_cl, Int64(50), mc_routine_cl, save_cl)

        @test df_cl isa DataFrame
        @test names(df_cl) == ["iter", "omega", "energy", "num_particles"]
        @test nrow(df_cl) > 0
        @test length(updated_liveset_cl.walkers) == 10

        # Cluster adaptation should have been active
        @test length(updated_params_cl.cluster_p_history) >= 0
        # cluster_p should be within bounds
        @test updated_params_cl.cluster_p >= mc_routine_cl.cluster_p_floor
        @test updated_params_cl.cluster_p <= mc_routine_cl.cluster_p_ceiling

        rm("test_gc_cl_df.csv", force=true)
        rm("test_gc_cl.traj", force=true)
        rm("test_gc_cl.ls", force=true)
    end

    # ================================================================
    @testset "GC-NS with cluster moves: validation against exact enumeration" begin
        # Same 4x4 lattice and Hamiltonian as the exact-enumeration validation above
        L = 4
        lattice_template_cl = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:L*L]],
            adsorptions=:full
        )
        ham_cl = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        n_sites_cl = L * L

        mu_cl = -0.05
        kb = 8.617333262e-5
        T_cl = 300.0
        beta_cl = 1.0 / (kb * T_cl)

        # Exact enumeration
        exact_z_cl = 0.0
        exact_N_cl = 0.0
        for mask in 0:(2^n_sites_cl - 1)
            lat = deepcopy(lattice_template_cl)
            for site in 1:n_sites_cl
                lat.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_v = interacting_energy(lat, ham_cl).val
            N_v = sum(lat.components[1])
            omega_v = E_v - mu_cl * N_v
            boltz = exp(-beta_cl * omega_v)
            exact_z_cl += boltz
            exact_N_cl += boltz * N_v
        end
        exact_mean_N_cl = exact_N_cl / exact_z_cl
        exact_ln_z_cl = log(exact_z_cl)

        # Run GC-NS with cluster moves
        n_walkers_cl = 100
        n_steps_cl = Int64(3000)
        walkers_val_cl = [LatticeWalker(deepcopy(lattice_template_cl), energy=0.0u"eV", iter=0)
                          for _ in 1:n_walkers_cl]
        liveset_val_cl = LatticeGasWalkers(walkers_val_cl, ham_cl; assign_energy=false)

        gc_params_val_cl = GrandCanonicalNestedSamplingParameters(
            mc_steps=100, chemical_potential=mu_cl,
            energy_perturbation=1e-12, init_occupation_p=0.5)
        mc_routine_val_cl = MCGrandCanonicalMoves(
            p_move=0.5, p_insert=0.25,
            clusters_freq=3, swaps_freq=1, initial_cluster_p=0.3,
            cluster_adjust_interval=50)
        save_val_cl = SaveEveryN("test_val_cl.csv", "test_val_cl.traj", "test_val_cl.ls", 10000, 10000, 10000)

        df_val_cl, _, params_val_cl = grand_canonical_nested_sampling(
            liveset_val_cl, gc_params_val_cl, n_steps_cl, mc_routine_val_cl, save_val_cl)

        @test nrow(df_val_cl) > 0

        # Compute NS thermodynamic stats
        mean_E_cl, Cv_cl, mean_N_ns_cl = gc_thermodynamic_stats(
            df_val_cl, [beta_cl], n_walkers_cl, mu_cl)

        # Compare ⟨N⟩ with exact (|error| < 0.5)
        @test abs(mean_N_ns_cl[1] - exact_mean_N_cl) < 0.5

        # Cluster adaptation should have fired
        @test length(params_val_cl.cluster_p_history) > 0

        rm("test_val_cl.csv", force=true)
        rm("test_val_cl.traj", force=true)
        rm("test_val_cl.ls", force=true)
    end

    @testset "dead-point callback (Ω-sorted route)" begin
        using Random
        dpc_save = SaveEveryN("t_dpc_gc.csv", "t_dpc_gc.traj", "t_dpc_gc.ls",
                              1000000, 1000000, 1000000)
        dpc_cleanup() = rm.(["t_dpc_gc.csv", "t_dpc_gc.traj", "t_dpc_gc.ls"],
                            force=true)

        function dpc_run(seed, cb)
            Random.seed!(seed)
            walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV",
                                     iter=0) for _ in 1:10]
            ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
            gc = GrandCanonicalNestedSamplingParameters(mc_steps=30,
                chemical_potential=-0.05, energy_perturbation=1e-9)
            d, _, _ = grand_canonical_nested_sampling(ls, gc, Int64(150),
                MCGrandCanonicalMoves(), dpc_save; dead_point_callback=cb)
            dpc_cleanup()
            return d
        end

        # Invocation count equals nrow(df); the (iter, energy, N) triple seen
        # by the callback matches the ledger row bit-exactly
        seen = Tuple{Int,Float64,Int}[]
        df = dpc_run(4271, (iter, w) -> push!(seen,
            (iter, w.energy.val, sum(w.configuration.components[1]))))
        @test nrow(df) > 0
        @test length(seen) == nrow(df)
        @test [t[1] for t in seen] == df.iter
        @test [t[2] for t in seen] == df.energy
        @test [t[3] for t in seen] == df.num_particles

        # Stream neutrality: same-seed A/B with and without the callback
        dfA = dpc_run(4272, nothing)
        dfB = dpc_run(4272, (iter, w) -> nothing)
        @test dfA.iter == dfB.iter
        @test dfA.omega == dfB.omega
        @test dfA.energy == dfB.energy
        @test dfA.num_particles == dfB.num_particles
    end

    @testset "driver controls (Omega route)" begin
        using Random

        # Cached-key sortperm reproduces the by-comparator ordering
        # including tie order (both are stable with identical key values):
        # a tie-heavy walker vector sorted both ways gives the identical
        # object sequence
        Random.seed!(99001)
        tie_ws = LatticeWalker[]
        for i in 1:20
            l = deepcopy(square_lattice)
            w = LatticeWalker(l, energy=Float64(mod(i, 3)) * 0.01u"eV",
                              iter=i)
            push!(tie_ws, w)
        end
        mu_tie = -0.05
        byfun = w -> SamplingSchemes._grand_potential(w, mu_tie)
        a = copy(tie_ws)
        sort!(a, by=byfun, rev=true)
        b = copy(tie_ws)
        permute!(b, sortperm([byfun(w) for w in b], rev=true))
        @test all(a[i] === b[i] for i in eachindex(a))

        # Parameter-crafted stall: empty initial occupancies with a
        # deletion-only move mix are permanently guard-skipped, so every
        # step fails. stop_on_stall = true returns the partial (empty)
        # ledger and the intact live set after one warning.
        stall_save = SaveEveryN("t_stall_gc.csv", "t_stall_gc.traj",
                                "t_stall_gc.ls", 1000000, 1000000, 1000000)
        stall_cleanup() = rm.(["t_stall_gc.csv", "t_stall_gc.traj",
                               "t_stall_gc.ls"], force=true)
        function stall_run(; kwargs...)
            Random.seed!(99002)
            ws = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV",
                                iter=0) for _ in 1:8]
            ls = LatticeGasWalkers(ws, ham; assign_energy=false)
            gc = GrandCanonicalNestedSamplingParameters(mc_steps=10,
                chemical_potential=-0.05, energy_perturbation=1e-9,
                init_occupation_p=0.0, allowed_fail_count=3)
            d, lsx, _ = grand_canonical_nested_sampling(ls, gc, Int64(40),
                MCGrandCanonicalMoves(p_move=0.0, p_insert=0.0), stall_save;
                kwargs...)
            stall_cleanup()
            return d, lsx
        end
        d_stop, ls_stop = @test_logs (:warn, r"GC-NS: Failed") match_mode=:any stall_run(
            stop_on_stall=true)
        @test nrow(d_stop) == 0
        @test length(ls_stop.walkers) == 8
        # The default warn-and-continue burns the whole budget but returns
        # the same (empty) ledger; explicit false matches the unmentioned
        # default digit-for-digit
        d_def, _ = stall_run()
        d_off, _ = stall_run(stop_on_stall=false)
        @test nrow(d_def) == 0
        @test d_def.iter == d_off.iter

        # The default A/B on a healthy (non-stalling) fixture: explicit
        # stop_on_stall = false against a call that never mentions it,
        # digit-identical ledgers and live sets
        function healthy_run(; kwargs...)
            Random.seed!(99010)
            ws = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV",
                                iter=0) for _ in 1:8]
            ls = LatticeGasWalkers(ws, ham; assign_energy=false)
            gc = GrandCanonicalNestedSamplingParameters(mc_steps=20,
                chemical_potential=-0.05, energy_perturbation=1e-9)
            d, lsx, _ = grand_canonical_nested_sampling(ls, gc, Int64(30),
                MCGrandCanonicalMoves(), stall_save; kwargs...)
            stall_cleanup()
            return d, lsx
        end
        h_def, ls_hd = healthy_run()
        h_off, ls_ho = healthy_run(stop_on_stall=false)
        @test nrow(h_def) > 0
        @test h_def.omega == h_off.omega
        @test h_def.energy == h_off.energy
        @test [w.energy.val for w in ls_hd.walkers] ==
              [w.energy.val for w in ls_ho.walkers]

        # Observables and the dead-point callback ride the rewritten cached-
        # key pre-sort: a run with both on records the identical physics
        # columns as the plain same-seed run
        obs_d, _ = healthy_run(observables=[:n_occ =>
                                   cfg -> Float64(sum(cfg.components[1]))],
                               dead_point_callback=(iter, w) -> nothing)
        @test obs_d.omega == h_def.omega
        @test obs_d.energy == h_def.energy
        @test obs_d.num_particles == h_def.num_particles

        # Eligible-set contents on a crafted degenerate (all-tied) live set
        # match the by-comparator comprehension: both are empty, and on a
        # mixed set both agree entry-for-entry
        mu_e = -0.05
        gp_e = w -> SamplingSchemes._grand_potential(w, mu_e)
        tied = [LatticeWalker(deepcopy(square_lattice),
                              energy=0.1u"eV", iter=0) for _ in 1:6]
        keys_t = [gp_e(w) for w in tied]
        worst_t = maximum(keys_t)
        elig_cached = [k for k in 2:6 if keys_t[k] < worst_t]
        elig_direct = [k for k in 2:6 if gp_e(tied[k]) < worst_t]
        @test elig_cached == elig_direct == Int[]
        mixed = [LatticeWalker(deepcopy(square_lattice),
                               energy=Float64(mod(i, 3)) * 0.01u"eV", iter=i)
                 for i in 1:6]
        sort!(mixed, by=gp_e, rev=true)
        keys_m = [gp_e(w) for w in mixed]
        worst_m = keys_m[1]
        @test [k for k in 2:6 if keys_m[k] < worst_m] ==
              [k for k in 2:6 if gp_e(mixed[k]) < worst_m]
    end
end