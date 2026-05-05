@testset "U-sorted ideal-gas-reference GCNS tests" begin

    # ================================================================
    # Shared lattice and Hamiltonian for all tests
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
    @testset "LatticeGCNSIdealGasReferenceParameters" begin
        p = LatticeGCNSIdealGasReferenceParameters()
        @test p isa SamplingSchemes.SamplingParameters
        @test p.mc_steps == 100
        @test p.z0 == 1.0
        @test p.energy_perturbation == 1e-12
        @test p.random_seed == 1234
        @test p.fail_count == 0
        @test p.allowed_fail_count == 10
        @test p.n_max == typemax(Int64)
        @test p.T0_ref === nothing

        # Cluster runtime state defaults
        @test p.cluster_p == 0.3
        @test p.cluster_accepted == 0.0
        @test p.cluster_total == 0.0
        @test p.cluster_p_history == Float64[]
        @test p.cluster_accept_history == Float64[]
        @test p.cluster_adjust_iterations == Int[]

        # Custom keyword args
        p2 = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=200, z0=2.5, n_max=Int64(8), T0_ref=300.0)
        @test p2.mc_steps == 200
        @test p2.z0 == 2.5
        @test p2.n_max == 8
        @test p2.T0_ref == 300.0

        # Mutability of runtime state
        p.fail_count = 7
        @test p.fail_count == 7
        p.cluster_p = 0.5
        @test p.cluster_p == 0.5

        # Validation
        @test_throws ArgumentError LatticeGCNSIdealGasReferenceParameters(z0=0.0)
        @test_throws ArgumentError LatticeGCNSIdealGasReferenceParameters(z0=-1.0)
        @test_throws ArgumentError LatticeGCNSIdealGasReferenceParameters(T0_ref=0.0)
        @test_throws ArgumentError LatticeGCNSIdealGasReferenceParameters(T0_ref=-100.0)
    end

    # ================================================================
    @testset "MCIdealGasReferenceMoves" begin
        mc = MCIdealGasReferenceMoves()
        @test mc isa MCRoutine
        @test mc.p_flip == 0.5
        @test mc.swaps_freq == 1
        @test mc.clusters_freq == 0
        @test mc.initial_cluster_p == 0.3
        @test mc.target_cluster_accept == 0.3
        @test mc.cluster_adjust_interval == 50
        @test mc.cluster_p_floor == 0.01
        @test mc.cluster_p_ceiling == 1.0

        # Custom args
        mc2 = MCIdealGasReferenceMoves(
            p_flip=0.4, swaps_freq=2, clusters_freq=3, initial_cluster_p=0.6)
        @test mc2.p_flip == 0.4
        @test mc2.swaps_freq == 2
        @test mc2.clusters_freq == 3
        @test mc2.initial_cluster_p == 0.6

        # Allow p_flip = 1 (no fixed-N moves required)
        mc3 = MCIdealGasReferenceMoves(p_flip=1.0, swaps_freq=0, clusters_freq=0)
        @test mc3.p_flip == 1.0

        # Validation
        @test_throws ArgumentError MCIdealGasReferenceMoves(p_flip=-0.1)
        @test_throws ArgumentError MCIdealGasReferenceMoves(p_flip=1.1)
        @test_throws ArgumentError MCIdealGasReferenceMoves(p_flip=0.5, swaps_freq=0, clusters_freq=0)
        @test_throws ArgumentError MCIdealGasReferenceMoves(swaps_freq=-1)
        @test_throws ArgumentError MCIdealGasReferenceMoves(clusters_freq=-1)
        @test_throws ArgumentError MCIdealGasReferenceMoves(initial_cluster_p=1.5)
        @test_throws ArgumentError MCIdealGasReferenceMoves(cluster_p_floor=0.5, cluster_p_ceiling=0.3)
    end

    # ================================================================
    @testset "_init_ideal_gas_reference_walkers!" begin
        n_sites = num_sites(square_lattice)  # 16

        # z0 = 1 → p0 = 0.5; mean N over many walkers should be close to n_sites/2
        walkers_1 = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                     for _ in 1:200]
        liveset_1 = LatticeGasWalkers(walkers_1, ham; assign_energy=false)
        params_1 = LatticeGCNSIdealGasReferenceParameters(z0=1.0, energy_perturbation=0.0)
        SamplingSchemes._init_ideal_gas_reference_walkers!(liveset_1, params_1)
        Ns_1 = [sum(w.configuration.components[1]) for w in liveset_1.walkers]
        @test mean(Ns_1) ≈ n_sites * 0.5 atol=1.0  # ~8 ± 1 over 200 walkers
        # Energies should be assigned (not zero unless lattice happens to be empty)
        @test all(w.energy isa typeof(0.0u"eV") for w in liveset_1.walkers)

        # z0 = 0.01 → p0 ≈ 0.0099; very low occupancy
        walkers_low = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                       for _ in 1:200]
        liveset_low = LatticeGasWalkers(walkers_low, ham; assign_energy=false)
        params_low = LatticeGCNSIdealGasReferenceParameters(z0=0.01)
        SamplingSchemes._init_ideal_gas_reference_walkers!(liveset_low, params_low)
        Ns_low = [sum(w.configuration.components[1]) for w in liveset_low.walkers]
        @test mean(Ns_low) < 1.0  # well under 1 expected (mean ≈ 16 * 0.0099 ≈ 0.16)

        # z0 = 100 → p0 ≈ 0.99; very high occupancy
        walkers_high = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                        for _ in 1:200]
        liveset_high = LatticeGasWalkers(walkers_high, ham; assign_energy=false)
        params_high = LatticeGCNSIdealGasReferenceParameters(z0=100.0)
        SamplingSchemes._init_ideal_gas_reference_walkers!(liveset_high, params_high)
        Ns_high = [sum(w.configuration.components[1]) for w in liveset_high.walkers]
        @test mean(Ns_high) > n_sites - 1.0  # well over n_sites - 1 expected (~15.84)

        # n_max enforcement
        walkers_cap = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                       for _ in 1:50]
        liveset_cap = LatticeGasWalkers(walkers_cap, ham; assign_energy=false)
        params_cap = LatticeGCNSIdealGasReferenceParameters(z0=10.0, n_max=Int64(5))
        SamplingSchemes._init_ideal_gas_reference_walkers!(liveset_cap, params_cap)
        for w in liveset_cap.walkers
            @test sum(w.configuration.components[1]) <= 5
        end

        # Energy is reassigned to match configuration after init
        walker_eng = LatticeWalker(deepcopy(square_lattice), energy=999.0u"eV", iter=0)
        liveset_eng = LatticeGasWalkers([walker_eng], ham; assign_energy=false)
        params_eng = LatticeGCNSIdealGasReferenceParameters(z0=1.0, energy_perturbation=0.0)
        SamplingSchemes._init_ideal_gas_reference_walkers!(liveset_eng, params_eng)
        recomputed = interacting_energy(liveset_eng.walkers[1].configuration, ham)
        @test liveset_eng.walkers[1].energy ≈ recomputed
    end

    # ================================================================
    @testset "MC_ideal_gas_reference_walk! basic" begin
        walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        assign_energy!(walker, ham)
        # Generous u_max so ceiling never binds (initial energy is negative
        # since on-site is -0.04, neighbors are negative; +1 is well above)
        u_max = walker.energy.val + 1.0

        accept, rate, updated, cl_acc, cl_tot = MC_ideal_gas_reference_walk!(
            100, walker, ham, u_max, 1.0;
            p_flip=0.5, swaps_freq=1, clusters_freq=0,
            energy_perturb=0.0)

        @test accept isa Bool
        @test 0.0 <= rate <= 1.0
        @test updated isa LatticeWalker{1}
        # No cluster moves when clusters_freq = 0
        @test cl_acc == 0
        @test cl_tot == 0

        # Stored energy matches recomputed energy of the stored configuration
        @test updated.energy ≈ interacting_energy(updated.configuration, ham)

        # Energy ceiling respected
        @test updated.energy.val < u_max

        # n_max respected
        @test sum(updated.configuration.components[1]) <= 16

        # Validation
        @test_throws ArgumentError MC_ideal_gas_reference_walk!(
            10, walker, ham, u_max, 1.0; p_flip=-0.1)
        @test_throws ArgumentError MC_ideal_gas_reference_walk!(
            10, walker, ham, u_max, 1.0; p_flip=1.1)
        @test_throws ArgumentError MC_ideal_gas_reference_walk!(
            10, walker, ham, u_max, 0.0)
        @test_throws ArgumentError MC_ideal_gas_reference_walk!(
            10, walker, ham, u_max, 1.0; p_flip=0.5, swaps_freq=0, clusters_freq=0)
    end

    # ================================================================
    @testset "MC_ideal_gas_reference_walk! with cluster moves" begin
        walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        assign_energy!(walker, ham)
        u_max = walker.energy.val + 1.0

        accept, rate, updated, cl_acc, cl_tot = MC_ideal_gas_reference_walk!(
            300, walker, ham, u_max, 1.0;
            p_flip=0.4, swaps_freq=1, clusters_freq=1, cluster_p=0.3,
            energy_perturb=0.0)

        @test accept isa Bool
        @test 0.0 <= rate <= 1.0
        @test updated isa LatticeWalker{1}
        @test cl_tot > 0
        @test 0 <= cl_acc <= cl_tot

        @test updated.energy ≈ interacting_energy(updated.configuration, ham)
    end

    # ================================================================
    @testset "MC_ideal_gas_reference_walk! n_max enforcement" begin
        walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
        # Start with a known small N
        walker.configuration.components[1] .= false
        walker.configuration.components[1][1] = true
        walker.configuration.components[1][2] = true
        assign_energy!(walker, ham)
        u_max = 1.0  # generous

        for _ in 1:5  # repeat a few times for statistical confidence
            walker.configuration.components[1] .= false
            walker.configuration.components[1][1] = true
            walker.configuration.components[1][2] = true
            assign_energy!(walker, ham)
            _, _, updated, _, _ = MC_ideal_gas_reference_walk!(
                500, walker, ham, u_max, 5.0;  # large z0 wants to grow
                p_flip=1.0, swaps_freq=0, clusters_freq=0,
                n_max=Int64(4), energy_perturb=0.0)
            @test sum(updated.configuration.components[1]) <= 4
        end
    end

    # ================================================================
    @testset "MC_ideal_gas_reference_walk! detailed balance — non-interacting prior" begin
        # With a zero Hamiltonian, U(σ) = 0 for all σ and the energy ceiling
        # never binds (any u_max > 0 is loose). The chain then samples the
        # ideal-gas reference prior π₀(σ) ∝ z₀^N(σ), which is a per-site
        # Bernoulli with p₀ = z₀/(1+z₀). Mean N should match n_sites · p₀.
        ham_zero = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        n_sites = num_sites(square_lattice)  # 16
        u_max = 1.0

        function run_mean_N(z0, n_walkers, mc_steps)
            Ns = Vector{Int}(undef, n_walkers)
            for k in 1:n_walkers
                w = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                # Random initial state at p=0.5 (uninformative)
                random_microstate!(w.configuration; p=0.5)
                assign_energy!(w, ham_zero)
                _, _, updated, _, _ = MC_ideal_gas_reference_walk!(
                    mc_steps, w, ham_zero, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                Ns[k] = sum(updated.configuration.components[1])
            end
            return mean(Ns)
        end

        # z0 = 1 → ⟨N⟩ = 8
        m1 = run_mean_N(1.0, 100, 200)
        @test m1 ≈ 8.0 atol=1.0

        # z0 = 2 → ⟨N⟩ = 16 · 2/3 ≈ 10.67
        m2 = run_mean_N(2.0, 100, 200)
        @test m2 ≈ 16.0 * 2.0/3.0 atol=1.0

        # z0 = 0.5 → ⟨N⟩ = 16 · 1/3 ≈ 5.33
        m05 = run_mean_N(0.5, 100, 200)
        @test m05 ≈ 16.0 / 3.0 atol=1.0
    end

    # ================================================================
    @testset "nested_sampling_step! for ideal-gas-reference GCNS" begin
        walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                   for _ in 1:5]
        liveset = LatticeGasWalkers(walkers, ham)

        params = LatticeGCNSIdealGasReferenceParameters(mc_steps=50, z0=1.5)
        mc_routine = MCIdealGasReferenceMoves()

        # Initialize walkers from prior
        SamplingSchemes._init_ideal_gas_reference_walkers!(liveset, params)

        e_type = typeof(walkers[1].energy)
        iter, u_worst, n_worst, updated_liveset, updated_params = nested_sampling_step!(
            liveset, params, mc_routine)

        @test iter isa Union{Missing,Int}
        @test u_worst isa Union{Missing,e_type}
        @test n_worst isa Union{Missing,Int}
        @test length(updated_liveset.walkers) == 5
        @test updated_params.fail_count >= 0
    end

    # ================================================================
    @testset "ideal_gas_reference_grand_canonical_nested_sampling loop" begin
        walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                   for _ in 1:10]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)

        params = LatticeGCNSIdealGasReferenceParameters(mc_steps=50, z0=1.5)
        mc_routine = MCIdealGasReferenceMoves()
        save_strategy = SaveEveryN(
            df_filename="test_igr_df.csv",
            wk_filename="test_igr.traj",
            ls_filename="test_igr.ls",
            n_traj=1000, n_snap=1000, n_info=1000)

        df, updated_liveset, updated_params = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset, params, Int64(20), mc_routine, save_strategy)

        @test df isa DataFrame
        @test names(df) == ["iter", "emax", "num_particles"]
        @test nrow(df) <= 20
        @test nrow(df) > 0
        @test eltype(df.iter) == Int
        @test eltype(df.emax) == Float64
        @test eltype(df.num_particles) == Int
        @test length(updated_liveset.walkers) == 10

        # U is monotonically non-increasing (each recorded U <= first; small slack
        # for energy-perturbation tie-breaker)
        if nrow(df) > 1
            for i in 2:nrow(df)
                @test df.emax[i] <= df.emax[1] + 1e-10
            end
        end

        # Particle counts in valid range
        @test all(0 .<= df.num_particles .<= num_sites(square_lattice))

        rm("test_igr_df.csv", force=true)
        rm("test_igr.traj", force=true)
        rm("test_igr.ls", force=true)
    end

    # ================================================================
    @testset "ideal-gas-reference GCNS with cluster moves" begin
        walkers_cl = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                      for _ in 1:10]
        liveset_cl = LatticeGasWalkers(walkers_cl, ham; assign_energy=false)

        params_cl = LatticeGCNSIdealGasReferenceParameters(mc_steps=50, z0=1.0)
        mc_routine_cl = MCIdealGasReferenceMoves(
            p_flip=0.5, swaps_freq=1, clusters_freq=1, initial_cluster_p=0.3,
            cluster_adjust_interval=20)
        save_cl = SaveEveryN(
            df_filename="test_igr_cl_df.csv",
            wk_filename="test_igr_cl.traj",
            ls_filename="test_igr_cl.ls",
            n_traj=10000, n_snap=10000, n_info=10000)

        df_cl, updated_liveset_cl, updated_params_cl = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset_cl, params_cl, Int64(50), mc_routine_cl, save_cl)

        @test df_cl isa DataFrame
        @test names(df_cl) == ["iter", "emax", "num_particles"]
        @test nrow(df_cl) > 0
        @test length(updated_liveset_cl.walkers) == 10

        # cluster_p stays within bounds
        @test mc_routine_cl.cluster_p_floor <=
              updated_params_cl.cluster_p <=
              mc_routine_cl.cluster_p_ceiling

        rm("test_igr_cl_df.csv", force=true)
        rm("test_igr_cl.traj", force=true)
        rm("test_igr_cl.ls", force=true)
    end

    # ================================================================
    @testset "ideal-gas-reference GCNS n_max enforcement" begin
        n_max_val = Int64(5)
        walkers_nm = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                      for _ in 1:10]
        liveset_nm = LatticeGasWalkers(walkers_nm, ham; assign_energy=false)

        params_nm = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=50, z0=10.0, n_max=n_max_val)  # large z0 wants to grow
        mc_routine_nm = MCIdealGasReferenceMoves()
        save_nm = SaveEveryN(
            df_filename="test_igr_nmax_df.csv",
            wk_filename="test_igr_nmax.traj",
            ls_filename="test_igr_nmax.ls",
            n_traj=10000, n_snap=10000, n_info=10000)

        df_nm, updated_liveset_nm, _ = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset_nm, params_nm, Int64(30), mc_routine_nm, save_nm)

        for w in updated_liveset_nm.walkers
            @test sum(w.configuration.components[1]) <= n_max_val
        end

        if nrow(df_nm) > 0
            @test all(df_nm.num_particles .<= n_max_val)
        end

        rm("test_igr_nmax_df.csv", force=true)
        rm("test_igr_nmax.traj", force=true)
        rm("test_igr_nmax.ls", force=true)
    end

    # ================================================================
    @testset "ideal_gas_reference_thermodynamic_stats basic" begin
        # Hand-crafted: 2 microstates, equal NS weights, z0 = 1
        # At z0 = 1, the formula reduces to the Ω-sorted GCNS reweighting,
        # so we can cross-check against gc_thermodynamic_stats with
        # Ω = E - μN.
        μ = -0.5
        ωi = [0.5, 0.5]
        Us = [0.0, -1.0]
        Ns = [0, 1]

        # At low β
        β = 1.0
        e1, n1, cv1, lnz1 = ideal_gas_reference_thermodynamic_stats(
            β, μ, ωi, Us, Ns; z0=1.0, kb=1.0)
        @test isfinite(e1) && isfinite(n1) && isfinite(cv1) && isfinite(lnz1)

        # Cross-check against gc_thermodynamic_stats at z0 = 1
        omegas = [Us[i] - μ * Ns[i] for i in 1:2]
        e1_ref, cv1_ref, n1_ref = gc_thermodynamic_stats(
            β, ωi, omegas, Us, Ns, μ; kb=1.0)
        @test e1 ≈ e1_ref
        @test n1 ≈ n1_ref
        @test cv1 ≈ cv1_ref

        # Validation
        @test_throws DimensionMismatch ideal_gas_reference_thermodynamic_stats(
            β, μ, ωi, [0.0, 1.0, 2.0], Ns; z0=1.0)
        @test_throws ArgumentError ideal_gas_reference_thermodynamic_stats(
            β, μ, ωi, Us, Ns; z0=0.0)
        @test_throws ArgumentError ideal_gas_reference_thermodynamic_stats(
            β, μ, ωi, Us, Ns; z0=-1.0)

        # Empty input → all NaN
        e_empty, n_empty, cv_empty, lnz_empty = ideal_gas_reference_thermodynamic_stats(
            β, μ, Float64[], Float64[], Int[]; z0=1.0)
        @test isnan(e_empty) && isnan(n_empty) && isnan(cv_empty) && isnan(lnz_empty)
    end

    # ================================================================
    @testset "ideal_gas_reference_thermodynamic_stats smoke at varied z0" begin
        # Basic smoke test that the scalar function runs cleanly at
        # different z0 values on synthetic data and returns finite outputs.
        # (Equivalence to gc_thermodynamic_stats at z0=1 is checked above.)
        μ = -0.1
        β = 1.0 / (8.617333262e-5 * 300.0)
        ωi = [0.3, 0.3, 0.2, 0.1, 0.05, 0.05]
        Us = [-1.0, -0.5, 0.0, 0.2, 0.3, 0.5]
        Ns = [4, 3, 2, 1, 1, 0]
        for z0 in (0.5, 1.0, 2.0, 5.0)
            e, n, cv_mu, lnz = ideal_gas_reference_thermodynamic_stats(
                β, μ, ωi, Us, Ns; z0=z0)
            @test isfinite(e) && isfinite(n) && isfinite(cv_mu) && isfinite(lnz)
        end
    end

    # ================================================================
    @testset "ideal_gas_reference_thermodynamic_stats DataFrame sweep" begin
        # Build a synthetic DataFrame mimicking GCNS output
        df = DataFrame(
            iter=[1, 2, 3, 4, 5],
            emax=[-1.0, -0.5, 0.0, 0.5, 1.0],
            num_particles=[5, 4, 3, 2, 1])

        μs = [-0.05, -0.02, 0.0]
        Ts = [200.0, 300.0, 500.0]

        result = ideal_gas_reference_thermodynamic_stats(
            df, μs, Ts, 10, 1.0)

        @test result isa NamedTuple
        @test haskey(result, :mean_E)
        @test haskey(result, :mean_N)
        @test haskey(result, :Cv_mu)
        @test haskey(result, :ln_Z_relative)
        @test size(result.mean_E) == (3, 3)
        @test size(result.mean_N) == (3, 3)
        @test size(result.Cv_mu) == (3, 3)
        @test size(result.ln_Z_relative) == (3, 3)
        @test all(isfinite, result.mean_E)
        @test all(isfinite, result.mean_N)
    end

    # ================================================================
    @testset "effective_sample_size_ideal_gas_reference basic" begin
        # When ωi are all equal and r_j = 1 (T = T0, z = z0), N_eff = N.
        ωi = fill(1.0, 7)
        Us = collect(0.0:6.0)
        Ns = collect(0:6)

        z0 = 2.0
        T = 300.0
        T0 = 300.0
        kb = 8.617333262e-5
        μ = (1.0 / (1.0 / (kb * T))) * log(z0)  # makes z(μ,T) == z0
        β = 1.0 / (kb * T)

        # μ chosen so z = z0 → r_j = 1 for all j
        Neff = effective_sample_size_ideal_gas_reference(
            β, μ, ωi, Us, Ns; z0=z0, T=T, T0=T0)
        @test Neff ≈ 7.0 rtol=1e-10

        # Reweighting far from sampling distribution should give Neff < N
        μ_far = μ + 5.0
        Neff_far = effective_sample_size_ideal_gas_reference(
            β, μ_far, ωi, Us, Ns; z0=z0, T=T, T0=T0)
        @test 0.0 < Neff_far < 7.0

        # Validation
        @test_throws DimensionMismatch effective_sample_size_ideal_gas_reference(
            β, μ, ωi, [0.0, 1.0], Ns; z0=z0, T=T)
        @test_throws ArgumentError effective_sample_size_ideal_gas_reference(
            β, μ, ωi, Us, Ns; z0=0.0, T=T)
        @test_throws ArgumentError effective_sample_size_ideal_gas_reference(
            β, μ, ωi, Us, Ns; z0=z0, T=-1.0)
        @test_throws ArgumentError effective_sample_size_ideal_gas_reference(
            β, μ, ωi, Us, Ns; z0=z0, T=T, T0=-1.0)

        # Empty input → NaN
        Neff_empty = effective_sample_size_ideal_gas_reference(
            β, μ, Float64[], Float64[], Int[]; z0=z0, T=T)
        @test isnan(Neff_empty)
    end

    # ================================================================
    @testset "effective_sample_size_ideal_gas_reference T0 default" begin
        # When T0 === nothing, default is T → no temperature reweighting,
        # only fugacity ratio. Verify by comparing T0=nothing to T0=T.
        ωi = [0.4, 0.3, 0.2, 0.1]
        Us = [-1.0, -0.5, 0.0, 0.5]
        Ns = [3, 2, 1, 0]
        z0 = 1.5
        T = 300.0
        kb = 8.617333262e-5
        μ = -0.05
        β = 1.0 / (kb * T)

        Neff_default = effective_sample_size_ideal_gas_reference(
            β, μ, ωi, Us, Ns; z0=z0, T=T)
        Neff_explicit = effective_sample_size_ideal_gas_reference(
            β, μ, ωi, Us, Ns; z0=z0, T=T, T0=T)
        @test Neff_default ≈ Neff_explicit rtol=1e-12
    end

    # ================================================================
    @testset "effective_sample_size_ideal_gas_reference DataFrame sweep" begin
        df = DataFrame(
            iter=[1, 2, 3, 4, 5],
            emax=[-1.0, -0.5, 0.0, 0.5, 1.0],
            num_particles=[5, 4, 3, 2, 1])

        μs = [-0.05, 0.0]
        Ts = [200.0, 300.0, 500.0]

        Neff_grid = effective_sample_size_ideal_gas_reference(
            df, μs, Ts, 10, 1.0)

        @test size(Neff_grid) == (2, 3)
        @test all(>(0), Neff_grid)
    end

    # ================================================================
    @testset "Tier 0: non-interacting lattice gas — analytic reference" begin
        # All Hamiltonian terms are zero, so U(σ) = 0 for every microstate
        # and the only nontrivial contribution to the partition function is
        # the prior. The closed-form reference is:
        #   Ξ(μ, T) = (1 + z(μ,T))^M
        #   ⟨N⟩(μ, T) = M · z(μ,T) / (1 + z(μ,T))
        #   ⟨E⟩ = 0,   C_{V,μ} = 0
        L = 4
        n_sites = L * L
        ham_zero = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        lattice_template = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:n_sites]],
            adsorptions=:full)

        z0 = 1.0
        n_walkers = 200
        n_ns_steps = Int64(2000)

        walkers_t0 = [LatticeWalker(deepcopy(lattice_template), energy=0.0u"eV", iter=0)
                      for _ in 1:n_walkers]
        liveset_t0 = LatticeGasWalkers(walkers_t0, ham_zero; assign_energy=false)

        # `allowed_fail_count` is set high here because non-interacting NS
        # naturally stalls once u_max approaches the perturbation lower bound:
        # the underlying U(σ) is constant, so once the perturbation noise floor
        # is exhausted there are no microstates left below u_max. This is
        # algorithmically expected and not a bug.
        params_t0 = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=100, z0=z0, energy_perturbation=1e-12,
            allowed_fail_count=10000)
        mc_routine_t0 = MCIdealGasReferenceMoves(
            p_flip=1.0, swaps_freq=0, clusters_freq=0)
        save_t0 = SaveEveryN(
            df_filename="test_tier0_df.csv",
            wk_filename="test_tier0.traj",
            ls_filename="test_tier0.ls",
            n_traj=10000, n_snap=10000, n_info=10000)

        df_t0, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset_t0, params_t0, n_ns_steps, mc_routine_t0, save_t0)

        @test nrow(df_t0) > 0

        # Sweep target μ at fixed T. The range is restricted to keep the
        # importance-weight variance manageable: at z₀=1 and μ=±0.025 (z≈0.38..2.6),
        # CV² of the importance weight (z/z₀)^N over Bin(M, 0.5) ≈ 17, so
        # the effective sample size after reweighting drops to ESS ≈ K/18.
        T = 300.0
        kb = 8.617333262e-5
        β = 1.0 / (kb * T)
        μs = [-0.025, -0.0125, 0.0, 0.0125, 0.025]
        result = ideal_gas_reference_thermodynamic_stats(df_t0, μs, [T], n_walkers, z0)

        # ⟨E⟩ = 0 exactly (no interactions, perturbation cancels to <1e-12)
        @test all(abs.(result.mean_E) .< 1e-6)
        # C_{V,μ} = 0 exactly (Var(U) = 0)
        @test all(abs.(result.Cv_mu) .< 1e-6)

        # ⟨N⟩ matches M · z/(1+z). The tolerance reflects the importance-weight
        # variance amplification: σ(⟨N⟩) ≈ sqrt(Var(N) · (1+CV²)/K) ≤ ~0.5
        # at the extremes for K=200. atol=2.0 ≈ 4σ.
        for (i, μ) in enumerate(μs)
            z = exp(β * μ)
            expected_N = n_sites * z / (1.0 + z)
            @test result.mean_N[i, 1] ≈ expected_N atol=2.0
        end

        # ln Ξ_absolute = M · ln(1 + z₀) + ln Z_NS_relative should match
        # M · ln(1 + z(μ,T)). Tolerance ~0.5 nat ≈ 5% of typical ln Ξ value
        # at the sweep extremes; for K=200 the NS estimator's SD on ln Ξ
        # is well under this.
        for (i, μ) in enumerate(μs)
            z = exp(β * μ)
            expected_lnΞ = n_sites * log(1.0 + z)
            actual_lnΞ = n_sites * log(1.0 + z0) + result.ln_Z_relative[i, 1]
            @test actual_lnΞ ≈ expected_lnΞ atol=0.5
        end

        rm("test_tier0_df.csv", force=true)
        rm("test_tier0.traj", force=true)
        rm("test_tier0.ls", force=true)
    end

    # ================================================================
    @testset "Tier 1: L=4 interacting lattice gas — exact enumeration" begin
        # 4×4 lattice with the standard FreeBird GCNS validation Hamiltonian:
        # ε_ads = -0.04 eV, ε_nn = -0.01 eV, ε_nnn = -0.0025 eV.
        # Exact enumeration over 2^16 = 65536 microstates is the ground truth.
        # NS run uses z₀ = 1 (μ₀ = k_B T · ln(z₀) = 0, half-filling prior).
        # Sweep target μ in a window centered on μ₀ = 0.
        L = 4
        n_sites = L * L
        ham_t1 = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        lattice_template_t1 = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:n_sites]],
            adsorptions=:full)

        T = 300.0
        kb = 8.617333262e-5
        β = 1.0 / (kb * T)
        z0 = 1.0
        μ_window = 0.04
        μs = [-μ_window, -μ_window/2, 0.0, μ_window/2, μ_window]

        # ----- Exact enumeration -----
        # Cache (E, N) per microstate so we can compute exact stats at every μ.
        E_states = Vector{Float64}(undef, 2^n_sites)
        N_states = Vector{Int}(undef, 2^n_sites)
        for mask in 0:(2^n_sites - 1)
            lat = deepcopy(lattice_template_t1)
            for site in 1:n_sites
                lat.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_states[mask + 1] = interacting_energy(lat, ham_t1).val
            N_states[mask + 1] = sum(lat.components[1])
        end

        function exact_stats(μ)
            log_w = .-β .* (E_states .- μ .* N_states)
            max_log = maximum(log_w)
            w = exp.(log_w .- max_log)
            z_rel = sum(w)
            ln_Ξ = max_log + log(z_rel)
            mean_E = sum(w .* E_states) / z_rel
            mean_E2 = sum(w .* E_states.^2) / z_rel
            mean_N = sum(w .* N_states) / z_rel
            mean_EN = sum(w .* E_states .* N_states) / z_rel
            var_E = mean_E2 - mean_E^2
            cov_EN = mean_EN - mean_E * mean_N
            Cv = kb * β^2 * (var_E - μ * cov_EN)
            return mean_E, mean_N, Cv, ln_Ξ
        end

        # ----- NS run -----
        n_walkers = 100
        n_ns_steps = Int64(3000)

        walkers_t1 = [LatticeWalker(deepcopy(lattice_template_t1), energy=0.0u"eV", iter=0)
                      for _ in 1:n_walkers]
        liveset_t1 = LatticeGasWalkers(walkers_t1, ham_t1; assign_energy=false)

        # `allowed_fail_count` raised so the NS warning does not fire when the
        # live set bunches near the ground state in late iterations — this is
        # algorithmically expected for a finite system, not a bug.
        params_t1 = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=100, z0=z0, energy_perturbation=1e-12,
            allowed_fail_count=10000)
        mc_routine_t1 = MCIdealGasReferenceMoves(
            p_flip=0.5, swaps_freq=1, clusters_freq=0)
        save_t1 = SaveEveryN(
            df_filename="test_tier1_df.csv",
            wk_filename="test_tier1.traj",
            ls_filename="test_tier1.ls",
            n_traj=10000, n_snap=10000, n_info=10000)

        df_t1, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset_t1, params_t1, n_ns_steps, mc_routine_t1, save_t1)

        @test nrow(df_t1) > 0

        result_t1 = ideal_gas_reference_thermodynamic_stats(df_t1, μs, [T], n_walkers, z0)

        # ----- Compare against exact enumeration -----
        # Generous tolerances on stochastic NS results, matching the existing
        # Ω-sorted GCNS validation conventions (rtol=0.3 on ⟨E⟩, ⟨N⟩).
        # The high-N_eff window is centered on μ₀ = 0; tolerances are looser
        # at the sweep extremes where importance-weight variance grows.
        for (i, μ) in enumerate(μs)
            ex_E, ex_N, ex_Cv, ex_lnΞ = exact_stats(μ)

            @test result_t1.mean_N[i, 1] ≈ ex_N rtol=0.3 atol=1.0
            @test result_t1.mean_E[i, 1] ≈ ex_E rtol=0.3 atol=0.05

            # Cv is hardest to converge; just check it has the right sign.
            if isfinite(result_t1.Cv_mu[i, 1]) && isfinite(ex_Cv) && ex_Cv > 0
                @test result_t1.Cv_mu[i, 1] > 0
            end

            # ln Ξ_absolute = M·ln(1 + z₀) + ln Z_NS_relative
            ns_lnΞ = n_sites * log(1.0 + z0) + result_t1.ln_Z_relative[i, 1]
            @test ns_lnΞ ≈ ex_lnΞ atol=1.0
        end

        rm("test_tier1_df.csv", force=true)
        rm("test_tier1.traj", force=true)
        rm("test_tier1.ls", force=true)
    end

    # ================================================================
    @testset "Tier 2 fast cross-check: U-sorted vs Ω-sorted at L=4" begin
        # Run both GCNS constructions on the same Hamiltonian and verify
        # they agree on ⟨N⟩ and ⟨E⟩ at a target μ inside the U-sorted run's
        # high-N_eff window. The comprehensive L=8 head-to-head benchmark
        # (with wall-time and (μ,T) coverage analysis) lives in
        # scripts/benchmark_gcns_ideal_gas_reference_vs_omega.jl.
        L = 4
        n_sites = L * L
        ham_t2 = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        lattice_template_t2 = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:n_sites]],
            adsorptions=:full)

        T = 300.0
        kb = 8.617333262e-5
        β = 1.0 / (kb * T)
        μ_target = 0.0  # at half-filling for z₀=1, well inside the U-sorted window

        n_walkers = 100
        n_steps = Int64(1500)

        # ----- Ω-sorted GCNS run (the existing construction) -----
        walkers_ω = [LatticeWalker(deepcopy(lattice_template_t2), energy=0.0u"eV", iter=0)
                     for _ in 1:n_walkers]
        liveset_ω = LatticeGasWalkers(walkers_ω, ham_t2; assign_energy=false)

        params_ω = GrandCanonicalNestedSamplingParameters(
            mc_steps=100, chemical_potential=μ_target,
            energy_perturbation=1e-12, allowed_fail_count=10000)
        routine_ω = MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25)
        save_ω = SaveEveryN(
            df_filename="test_tier2_omega_df.csv",
            wk_filename="test_tier2_omega.traj",
            ls_filename="test_tier2_omega.ls",
            n_traj=10000, n_snap=10000, n_info=10000)

        df_ω, _, _ = grand_canonical_nested_sampling(
            liveset_ω, params_ω, n_steps, routine_ω, save_ω)

        ω_mean_E, ω_Cv, ω_mean_N = gc_thermodynamic_stats(
            df_ω, [β], n_walkers, μ_target)

        # ----- U-sorted, ideal-gas-reference GCNS run -----
        walkers_u = [LatticeWalker(deepcopy(lattice_template_t2), energy=0.0u"eV", iter=0)
                     for _ in 1:n_walkers]
        liveset_u = LatticeGasWalkers(walkers_u, ham_t2; assign_energy=false)

        z0 = 1.0  # → μ₀ = 0, matching μ_target
        params_u = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=100, z0=z0,
            energy_perturbation=1e-12, allowed_fail_count=10000)
        routine_u = MCIdealGasReferenceMoves(
            p_flip=0.5, swaps_freq=1, clusters_freq=0)
        save_u = SaveEveryN(
            df_filename="test_tier2_u_df.csv",
            wk_filename="test_tier2_u.traj",
            ls_filename="test_tier2_u.ls",
            n_traj=10000, n_snap=10000, n_info=10000)

        df_u, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset_u, params_u, n_steps, routine_u, save_u)

        result_u = ideal_gas_reference_thermodynamic_stats(
            df_u, [μ_target], [T], n_walkers, z0)
        u_mean_E = result_u.mean_E[1, 1]
        u_mean_N = result_u.mean_N[1, 1]
        u_Cv = result_u.Cv_mu[1, 1]

        # ----- Cross-check -----
        # Both estimators target the same physical ⟨E⟩, ⟨N⟩, Cv. Differences
        # within statistical error are expected; tolerances accommodate the
        # combined error of both runs.
        @test u_mean_N ≈ ω_mean_N[1] rtol=0.3 atol=1.5
        @test u_mean_E ≈ ω_mean_E[1] rtol=0.3 atol=0.1

        # Both Cv estimates should be positive (a system at half-filling with
        # attractive interactions has nonzero specific heat).
        if isfinite(ω_Cv[1]) && isfinite(u_Cv)
            @test ω_Cv[1] > 0
            @test u_Cv > 0
        end

        rm("test_tier2_omega_df.csv", force=true)
        rm("test_tier2_omega.traj", force=true)
        rm("test_tier2_omega.ls", force=true)
        rm("test_tier2_u_df.csv", force=true)
        rm("test_tier2_u.traj", force=true)
        rm("test_tier2_u.ls", force=true)
    end

end
