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
end
