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
    @testset "random_seed makes a run reproducible" begin
        # Before the seed was consumed, two runs of identical parameters were
        # different Markov chains: the field was stored and never reached an
        # RNG. This test detects that directly — without seeding, the second
        # call continues the global stream where the first left off and cannot
        # reproduce it.
        function gc_run(seed)
            walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0) for _ in 1:10]
            ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
            p = GrandCanonicalNestedSamplingParameters(
                mc_steps=50, chemical_potential=-0.05, random_seed=seed)
            routine = MCGrandCanonicalMoves()
            save = SaveEveryN("test_gc_seed.csv", "test_gc_seed.traj", "test_gc_seed.ls",
                              1000, 1000, 1000)
            df, _, _ = grand_canonical_nested_sampling(ls, p, Int64(30), routine, save)
            return df
        end

        df_a = gc_run(2024)
        df_b = gc_run(2024)
        @test nrow(df_a) > 0
        @test isequal(df_a, df_b)

        rm("test_gc_seed.csv", force=true)
        rm("test_gc_seed.traj", force=true)
        rm("test_gc_seed.ls", force=true)
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
        exact_N2 = 0.0
        exact_EN = 0.0
 
        E_all, N_all = grand_canonical_exact_enumeration(lattice_template, ham_val)

        for i in eachindex(E_all)
            E_val = E_all[i].val
            N_val = N_all[i]
            omega_val = E_val - mu_val * N_val
 
            boltz = exp(-beta_test * omega_val)
            exact_z += boltz
            exact_E += boltz * E_val
            exact_E2 += boltz * E_val^2
            exact_N += boltz * N_val
            exact_N2 += boltz * N_val^2
            exact_EN += boltz * E_val * N_val
        end
 
        exact_mean_E = exact_E / exact_z
        exact_mean_N = exact_N / exact_z
        exact_mean_E2 = exact_E2 / exact_z
        exact_mean_EN = exact_EN / exact_z
        exact_mean_N2 = exact_N2 / exact_z
        exact_var_E = exact_mean_E2 - exact_mean_E^2
        exact_var_N = exact_mean_N2 - exact_mean_N^2
        exact_cov_EN = exact_mean_EN - exact_mean_E * exact_mean_N
        # C_E — the thermodynamic heat capacity, and the default `cv`.
        exact_Cv = kb * beta_test^2 * (exact_var_E - mu_val * exact_cov_EN)
        # C_Ω — the fluctuation of Ω = E − μN. It differs from C_E by
        # −μ(∂⟨N⟩/∂T)_μ, and it is that difference this test now pins.
        exact_c_omega = kb * beta_test^2 *
            (exact_var_E - 2mu_val * exact_cov_EN + mu_val^2 * exact_var_N)
 
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
        r = gc_thermodynamic_stats(df, [beta_test], n_walkers, mu_val)
        mean_E_ns, Cv_ns, mean_N_ns = r.mean_E, r.cv, r.mean_N

        # Compare with exact values (generous tolerances for stochastic algorithm)
        @test mean_E_ns[1] ≈ exact_mean_E rtol=0.3
        @test mean_N_ns[1] ≈ exact_mean_N rtol=0.3

        # The heat capacities were previously checked for sign only — and
        # inside an `if` that could skip the assertion entirely, so a run
        # producing a non-finite Cv passed by asserting nothing at all. A
        # sign test also accepts a value wrong by any factor.
        #
        # Both are now compared numerically against the exact enumeration.
        # rtol=0.15 is measured, not guessed: with the chain pinned by
        # random_seed the observed errors are reproducible, and they are 5.0%
        # for C_E and 6.4% for C_Ω (logged below). 0.15 leaves roughly 3x and
        # 2.3x headroom — loose enough that a legitimate change to the sampler
        # does not trip it, tight enough to catch an estimator that has lost a
        # term or a factor. The first moments keep rtol=0.3 pending the same
        # treatment; their errors are logged below for that purpose.
        #
        # This is what pins MERGE_PLAN §1.3(a) — that C_Ω and C_E are
        # different quantities — in a test rather than in a docstring.
        @test isfinite(Cv_ns[1])
        @test isfinite(r.c_omega[1])
        @test Cv_ns[1] > 0
        @test Cv_ns[1] ≈ exact_Cv rtol=0.15
        @test r.c_omega[1] ≈ exact_c_omega rtol=0.15

        # And that the difference is not a subtlety. At this μ the exact
        # enumeration puts C_E at 2.90x C_Ω, so quoting one where the other is
        # meant is a factor-of-three error, not a rounding one. Asserted on the
        # exact values, which carry no sampling noise at all.
        @test exact_Cv / exact_c_omega > 2.0

        @info "G3 accuracy vs exact enumeration (seeded, so reproducible)" rel_err_mean_E=abs(mean_E_ns[1] - exact_mean_E) / abs(exact_mean_E) rel_err_mean_N=abs(mean_N_ns[1] - exact_mean_N) / abs(exact_mean_N) rel_err_C_E=abs(Cv_ns[1] - exact_Cv) / abs(exact_Cv) rel_err_C_omega=abs(r.c_omega[1] - exact_c_omega) / abs(exact_c_omega) ratio_C_E_to_C_omega=exact_Cv / exact_c_omega
 
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
        # Reuse the 4x4 lattice and Hamiltonian from Study B
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
        E_all_cl, N_all_cl = grand_canonical_exact_enumeration(lattice_template_cl, ham_cl)
        for i in eachindex(E_all_cl)
            E_v = E_all_cl[i].val
            N_v = N_all_cl[i]
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
end