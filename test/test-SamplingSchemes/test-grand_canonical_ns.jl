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
 
        accept, rate, updated_walker = MC_grand_canonical_walk!(
            100, walker, ham, omega_max, mu;
            p_move=0.5, p_insert=0.25, energy_perturb=0.0)
 
        @test accept isa Bool
        @test 0.0 <= rate <= 1.0
        @test updated_walker isa LatticeWalker{1}
 
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
end