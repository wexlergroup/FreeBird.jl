using Random
using Combinatorics: combinations

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
    # Sampler bug-localization tests at z₀ ≠ 1
    # (campaign artifacts; Step A's z₀=2 regression target is expected to
    # fail until the U-sorted ideal-gas-reference NS sampler bug is fixed.)
    # ================================================================

    """
    Exact grand-canonical ⟨N⟩(μ, T) on an L=4 lattice. Enumerates all 2^16
    microstates by bitmask iteration; mirrors the Ω-sorted exact-enum
    reference at `test-grand_canonical_ns.jl:306`.
    """
    function _exact_gc_mean_N(μ::Float64, T::Float64,
                              h::ClassicalHamiltonian,
                              lattice_template)
        n_sites = num_sites(lattice_template)
        β = 1.0 / (8.617333262e-5 * T)
        Z = 0.0
        sum_N = 0.0
        for mask in 0:(2^n_sites - 1)
            lattice = deepcopy(lattice_template)
            for site in 1:n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E = interacting_energy(lattice, h).val
            N = sum(lattice.components[1])
            w = exp(-β * (E - μ * N))
            Z += w
            sum_N += w * N
        end
        return sum_N / Z
    end

    """
    Exact ⟨N⟩ under the U-ceiling–restricted prior
    π₀(σ) ∝ z₀^N(σ) · 𝟙[U(σ) < u_max]. This is the equilibrium that
    `MC_ideal_gas_reference_walk!` is supposed to preserve at fixed u_max.
    """
    function _exact_constrained_mean_N(z0::Float64, u_max::Float64,
                                       h::ClassicalHamiltonian,
                                       lattice_template)
        n_sites = num_sites(lattice_template)
        Z = 0.0
        sum_N = 0.0
        for mask in 0:(2^n_sites - 1)
            lattice = deepcopy(lattice_template)
            for site in 1:n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E = interacting_energy(lattice, h).val
            if E < u_max
                N = sum(lattice.components[1])
                w = z0 ^ N
                Z += w
                sum_N += w * N
            end
        end
        return sum_N / Z
    end

    """
    Exact N-distribution P(N) under π₀(σ) ∝ z₀^N(σ) · 𝟙[U(σ) < u_max].
    Returns Vector{Float64} of length n_sites+1 (index N+1 → P(N)).
    Used by deep-compression tests for distribution-level (not just ⟨N⟩)
    diagnostic.
    """
    function _exact_constrained_N_distribution(z0::Float64, u_max::Float64,
                                                h::ClassicalHamiltonian,
                                                lattice_template)
        n_sites = num_sites(lattice_template)
        weights_by_N = zeros(Float64, n_sites + 1)
        for mask in 0:(2^n_sites - 1)
            lattice = deepcopy(lattice_template)
            for site in 1:n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E = interacting_energy(lattice, h).val
            if E < u_max
                N = sum(lattice.components[1])
                weights_by_N[N + 1] += z0 ^ N
            end
        end
        Z = sum(weights_by_N)
        return weights_by_N ./ Z
    end

    """
    Exact L=4 ground state U for the standard Hamiltonian. Computed once
    by enumeration; reused by deep-compression tests.
    """
    function _l4_ground_state_U(h::ClassicalHamiltonian, lattice_template)
        n_sites = num_sites(lattice_template)
        U_min = Inf
        for mask in 0:(2^n_sites - 1)
            lattice = deepcopy(lattice_template)
            for site in 1:n_sites
                lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E = interacting_energy(lattice, h).val
            if E < U_min
                U_min = E
            end
        end
        return U_min
    end

    """
    Exact ⟨N⟩ AND P(N) under the U-ceiling-restricted prior at L=6,
    enumerated over N ∈ [N_lo, M] only. At L=6 the full 2^36 ≈ 6.87e10
    enumeration is infeasible, but for u_max near the ground state, only
    high-N configs satisfy U < u_max — N_lo=30 covers all such configs at
    u_max ≈ ground+0.2. Returns (mean_N, PN_vector_of_length_M+1).
    """
    function _exact_constrained_distribution_l6(z0::Float64, u_max::Float64,
                                                 h::ClassicalHamiltonian,
                                                 lattice_template;
                                                 N_lo::Int=30)
        n_sites = num_sites(lattice_template)
        weights_by_N = zeros(Float64, n_sites + 1)
        sum_N = 0.0
        Z = 0.0
        # Validation: at N_lo the maximum-attractiveness arrangement gives
        # the minimum U; if min_U(N_lo) > u_max, N_lo configs don't reach
        # the constrained support and N_lo could be raised. We do NOT
        # check the converse (min_U(N_lo - 1) > u_max) — the caller must
        # pick N_lo low enough.
        for N in N_lo:n_sites
            z0_pow_N = z0 ^ N
            for combo in combinations(1:n_sites, N)
                lat = deepcopy(lattice_template)
                @inbounds for s in 1:n_sites
                    lat.components[1][s] = false
                end
                for s in combo
                    lat.components[1][s] = true
                end
                E = interacting_energy(lat, h).val
                if E < u_max
                    weights_by_N[N + 1] += z0_pow_N
                    sum_N += z0_pow_N * N
                    Z += z0_pow_N
                end
            end
        end
        return sum_N / Z, weights_by_N ./ Z
    end

    """
    L=6 ground state U for the standard Hamiltonian. Single config (all
    sites occupied); no enumeration needed.
    """
    function _l6_ground_state_U(h::ClassicalHamiltonian)
        L6 = MLattice{1,SquareLattice}(
            lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(6, 6, 1), periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[true for _ in 1:36]],
            adsorptions=:full)
        return interacting_energy(L6, h).val
    end

    # ================================================================
    @testset "Sampler at z₀=1: L=4 exact-enum agreement (control)" begin
        # Control test: the U-sorted ideal-gas-reference NS at z₀=1 should
        # reproduce exact-enum ⟨N⟩(μ, T) to within NS sampling noise.
        # Mirrors the Ω-sorted exact-enum pattern at
        # test-grand_canonical_ns.jl:306. Expected to pass.
        z0 = 1.0
        n_walkers = 100
        n_steps = Int64(10000)
        mc_steps = 100

        μ_grid = [-0.04, 0.0, 0.02]
        T_grid = [200.0, 300.0]

        exact_N = [_exact_gc_mean_N(μ, T, ham, square_lattice)
                   for μ in μ_grid, T in T_grid]

        Random.seed!(4001)
        walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                   for _ in 1:n_walkers]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=mc_steps, z0=z0, energy_perturbation=1e-12,
            allowed_fail_count=100000)
        routine = MCIdealGasReferenceMoves(p_flip=0.5, swaps_freq=1, clusters_freq=0)
        save = SaveEveryN(df_filename="test_sampler_z0_1.csv",
                          wk_filename="test_sampler_z0_1.traj",
                          ls_filename="test_sampler_z0_1.ls",
                          n_traj=10*n_steps, n_snap=10*n_steps, n_info=10*n_steps)
        df, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset, params, n_steps, routine, save)
        rm("test_sampler_z0_1.csv", force=true)
        rm("test_sampler_z0_1.traj", force=true)
        rm("test_sampler_z0_1.ls", force=true)

        stats = ideal_gas_reference_thermodynamic_stats(
            df, μ_grid, T_grid, n_walkers, z0)

        @test all(abs.(stats.mean_N .- exact_N) .≤ 0.2)
    end

    # ================================================================
    @testset "Sampler at z₀=2: L=4 exact-enum agreement (regression target)" begin
        # The L=4 reproduction of the L=8 single-run-z₀=2 bias against
        # exact enum, with the same (μ, T) coverage as the d2b campaign.
        # Expected to FAIL until the sampler bug is fixed.
        # Note: this test failing is the campaign artifact. Do not relax
        # the tolerance to make it pass — see the Diagnostic Failure
        # Protocol in the campaign prompt.
        z0 = 2.0
        n_walkers = 100
        n_steps = Int64(10000)
        mc_steps = 100

        μ_grid = [-0.04, 0.0, 0.02]
        T_grid = [200.0, 300.0]

        exact_N = [_exact_gc_mean_N(μ, T, ham, square_lattice)
                   for μ in μ_grid, T in T_grid]

        Random.seed!(4002)
        walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                   for _ in 1:n_walkers]
        liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = LatticeGCNSIdealGasReferenceParameters(
            mc_steps=mc_steps, z0=z0, energy_perturbation=1e-12,
            allowed_fail_count=100000)
        routine = MCIdealGasReferenceMoves(p_flip=0.5, swaps_freq=1, clusters_freq=0)
        save = SaveEveryN(df_filename="test_sampler_z0_2.csv",
                          wk_filename="test_sampler_z0_2.traj",
                          ls_filename="test_sampler_z0_2.ls",
                          n_traj=10*n_steps, n_snap=10*n_steps, n_info=10*n_steps)
        df, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
            liveset, params, n_steps, routine, save)
        rm("test_sampler_z0_2.csv", force=true)
        rm("test_sampler_z0_2.traj", force=true)
        rm("test_sampler_z0_2.ls", force=true)

        stats = ideal_gas_reference_thermodynamic_stats(
            df, μ_grid, T_grid, n_walkers, z0)

        @test all(abs.(stats.mean_N .- exact_N) .≤ 0.2)
    end

    # ================================================================
    @testset "B1: _init_ideal_gas_reference_walkers! preserves p₀ at z₀≠1" begin
        # Tight quantitative check at z₀ ∈ {0.5, 1.0, 2.0}: empirical
        # per-site occupancy must match p₀ = z₀/(1+z₀) to a 5σ binomial
        # envelope over 1000 walkers × 16 sites = 16,000 sites. The
        # existing line-95 testset only one-sidedly bounds occupancy at
        # extreme z₀; this closes that gap at the values that matter for
        # the campaign.
        n_walkers = 1000
        n_sites = num_sites(square_lattice)

        for z0 in (0.5, 1.0, 2.0)
            Random.seed!(8000 + round(Int, 1000 * z0))
            walkers = [LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
                       for _ in 1:n_walkers]
            liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
            params = LatticeGCNSIdealGasReferenceParameters(
                z0=z0, energy_perturbation=0.0)
            SamplingSchemes._init_ideal_gas_reference_walkers!(liveset, params)

            total_occ = sum(sum(w.configuration.components[1])
                            for w in liveset.walkers)
            p_emp = total_occ / (n_walkers * n_sites)
            p_true = z0 / (1.0 + z0)
            se = sqrt(p_true * (1 - p_true) / (n_walkers * n_sites))
            @test abs(p_emp - p_true) ≤ 5 * se
        end
    end

    # ================================================================
    @testset "B2: MC kernel preserves π₀×𝟙[U<u_max] at z₀=1 (control)" begin
        # Control: at z₀=1 the Metropolis ratio collapses to 1, so the
        # gate-ordering question doesn't apply — total rejection equals
        # U-ceiling rejection. Confirms the harness samples the
        # constrained prior correctly when no z₀-correction is involved.
        #
        # No methodological-gate assertion here: at z₀=1 the typical
        # equilibrium N is lower than at z₀=2 (8 vs 10.7 on M=16), so the
        # same u_max=-0.20 produces a much smaller U-ceiling rejection
        # rate (~2% vs ~33%, measured during the localization session).
        # The gate is designed for the z₀=2 case where it both binds and
        # could expose the gate-ordering hypothesis.
        z0 = 1.0
        u_max = -0.20

        exact_mean_N = _exact_constrained_mean_N(z0, u_max, ham, square_lattice)

        n_chains = 50
        n_chunks = 50
        steps_per_chunk = 100
        Ns = Int[]

        Random.seed!(8001)
        p0 = z0 / (1.0 + z0)
        for _ in 1:n_chains
            walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
            random_microstate!(walker.configuration; p=p0)
            assign_energy!(walker, ham)
            attempts = 0
            while walker.energy.val >= u_max && attempts < 10000
                random_microstate!(walker.configuration; p=p0)
                assign_energy!(walker, ham)
                attempts += 1
            end
            @assert walker.energy.val < u_max  "Initial walker exceeds u_max — try a less restrictive ceiling"

            for _ in 1:n_chunks
                _, _, walker, _, _ = MC_ideal_gas_reference_walk!(
                    steps_per_chunk, walker, ham, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                push!(Ns, sum(walker.configuration.components[1]))
            end
        end

        @test abs(mean(Ns) - exact_mean_N) ≤ 0.4
    end

    # ================================================================
    @testset "B2: MC kernel preserves π₀×𝟙[U<u_max] at z₀=2" begin
        # Gate-ordering hypothesis test: in MC_ideal_gas_reference_walk!,
        # the U-ceiling gate (line 1029 of random_walks.jl) precedes the
        # Metropolis gate (line 1033). Detailed balance against
        # π₀(σ) · 𝟙[U(σ) < u_max] with the z₀-Bernoulli prior requires the
        # U-rejection to be symmetric under σ ↔ σ′ — which is true for
        # ΔN=0 moves but NOT obviously so for ΔN=±1 flips when the post-
        # flip configuration's U-passability differs from the pre-flip
        # configuration's.
        #
        # The line-230 detailed-balance test passes at z₀=2 only because
        # u_max=1.0 never binds against H_zero. This test exercises the
        # same kernel with the standard interacting Hamiltonian and
        # u_max=-0.20 eV, chosen to bind non-trivially: the U-ceiling
        # rejection rate at z₀=1 (control above) provides a calibrated
        # measurement of the gate-exercise level.
        z0 = 2.0
        u_max = -0.20

        exact_mean_N = _exact_constrained_mean_N(z0, u_max, ham, square_lattice)

        n_chains = 50
        n_chunks = 50
        steps_per_chunk = 100
        Ns = Int[]
        chain_reject_rates = Float64[]

        Random.seed!(8002)
        p0 = z0 / (1.0 + z0)
        for _ in 1:n_chains
            walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
            random_microstate!(walker.configuration; p=p0)
            assign_energy!(walker, ham)
            attempts = 0
            while walker.energy.val >= u_max && attempts < 10000
                random_microstate!(walker.configuration; p=p0)
                assign_energy!(walker, ham)
                attempts += 1
            end
            @assert walker.energy.val < u_max  "Initial walker exceeds u_max — try a less restrictive ceiling"

            chain_accepts = 0
            chain_total = 0
            for _ in 1:n_chunks
                _, accept_rate, walker, _, _ = MC_ideal_gas_reference_walk!(
                    steps_per_chunk, walker, ham, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                chain_accepts += round(Int, accept_rate * steps_per_chunk)
                chain_total += steps_per_chunk
                push!(Ns, sum(walker.configuration.components[1]))
            end
            push!(chain_reject_rates, 1.0 - chain_accepts / chain_total)
        end

        # Methodological gate: at z₀=2, total rejection includes both
        # U-ceiling and Metropolis. ≥5% per chain is the weakest possible
        # condition guaranteeing the kernel is doing real work — a tighter
        # condition on U-ceiling rejection specifically requires
        # instrumentation (see campaign report for the per-session
        # measurement of the U-ceiling component).
        @test minimum(chain_reject_rates) ≥ 0.05

        # Equilibrium check — the localizing assertion. If gate-ordering
        # is correct, empirical ⟨N⟩ matches the analytical restricted
        # prior to the MCMC SE (≤ 0.2 at this sample budget). If wrong,
        # the bias from L=8 (~2.5 particles in 64 sites) scales to ~0.6
        # particles at L=4 — well outside this tolerance.
        @test abs(mean(Ns) - exact_mean_N) ≤ 0.4
    end

    # ================================================================
    @testset "B2 deep-compression: MC kernel at z₀=1, u_max ≈ ground+0.1 (control)" begin
        # Control for the deep-compression test below. At z₀=1 the
        # Metropolis ratio collapses to 1, so the gate-ordering question
        # doesn't apply and the kernel reduces to U-ceiling-only filtering.
        # If z₀=1 here matches the analytical equilibrium under the same
        # u_max but the z₀=2 case below does not, the deviation is
        # specific to the z₀-Metropolis correction interacting with the
        # U-ceiling at near-ground-state regimes.
        z0 = 1.0
        ground_U = _l4_ground_state_U(ham, square_lattice)
        u_max = ground_U + 0.1   # tight ceiling, ~0.1 eV above ground state

        exact_mean_N = _exact_constrained_mean_N(z0, u_max, ham, square_lattice)
        exact_PN     = _exact_constrained_N_distribution(z0, u_max, ham, square_lattice)

        # Step-by-step observation so we can classify each accepted move
        # as insert (ΔN=+1) or delete (ΔN=-1). Per-step proposal type for
        # rejected moves is estimated from the per-step (M-N)/M and N/M
        # probabilities (uniform site selection).
        n_chains = 50
        n_steps_per_chain = 5000   # 250,000 total
        Ns = Int[]
        n_sites = num_sites(square_lattice)
        p0 = z0 / (1.0 + z0)

        ins_acc_total = 0
        del_acc_total = 0
        ins_total_estimated = 0.0
        del_total_estimated = 0.0

        Random.seed!(8101)
        for _ in 1:n_chains
            walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
            random_microstate!(walker.configuration; p=p0)
            assign_energy!(walker, ham)
            attempts = 0
            while walker.energy.val >= u_max && attempts < 50000
                random_microstate!(walker.configuration; p=p0)
                assign_energy!(walker, ham)
                attempts += 1
            end
            @assert walker.energy.val < u_max  "Initial walker exceeds u_max — adjust ceiling"

            for _ in 1:n_steps_per_chain
                N_before = sum(walker.configuration.components[1])
                ins_total_estimated += (n_sites - N_before) / n_sites
                del_total_estimated += N_before / n_sites
                _, _, walker, _, _ = MC_ideal_gas_reference_walk!(
                    1, walker, ham, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                N_after = sum(walker.configuration.components[1])
                if N_after - N_before == 1
                    ins_acc_total += 1
                elseif N_after - N_before == -1
                    del_acc_total += 1
                end
                push!(Ns, N_after)
            end
        end

        # Empirical N-distribution
        emp_PN = zeros(Float64, n_sites + 1)
        for N in Ns
            emp_PN[N + 1] += 1
        end
        emp_PN ./= length(Ns)
        tvd = 0.5 * sum(abs.(emp_PN .- exact_PN))

        # Per-class acceptance rates
        ins_acc_rate = ins_total_estimated > 0 ? ins_acc_total / ins_total_estimated : NaN
        del_acc_rate = del_total_estimated > 0 ? del_acc_total / del_total_estimated : NaN

        # Diagnostics (printed for inspection; not asserted)
        @info "B2 deep-compression z₀=1 control" u_max=u_max ground_U=ground_U exact_mean_N=exact_mean_N empirical_mean_N=mean(Ns) tvd=tvd ins_acc_rate=ins_acc_rate del_acc_rate=del_acc_rate ins_acc=ins_acc_total del_acc=del_acc_total ins_est=ins_total_estimated del_est=del_total_estimated

        # Equilibrium assertions — control should pass cleanly
        @test abs(mean(Ns) - exact_mean_N) ≤ 0.05
        @test tvd ≤ 0.05
    end

    # ================================================================
    @testset "B2 deep-compression: MC kernel at z₀=2, u_max ≈ ground+0.1" begin
        # Test the gate-ordering / insert-delete asymmetry hypothesis at the
        # L=4 analog of L=6's transitional regime, where the L=6 NS bias
        # accumulates (xi 13-18 contributing +0.40 to (μ,T)-reweighted ⟨N⟩).
        #
        # Setup: u_max ≈ ground_state + 0.1 eV. The constrained prior
        # π₀(σ) · 𝟙[U(σ) < u_max] is concentrated on configurations within
        # ~0.1 eV of the L=4 ground state — typically N ∈ {15, 16}, the
        # near-full-occupancy manifold where the U-ceiling binds tightly
        # against ΔN=+1 (insert) vs ΔN=-1 (delete) proposals asymmetrically.
        #
        # Hypothesis: at z₀=2, the kernel's insert/delete acceptance ratio
        # is mismatched against the constrained prior, biasing the
        # equilibrium ⟨N⟩ away from the analytical reference.
        #
        # If empirical ⟨N⟩ deviates from exact in the SAME direction as
        # L=6's transitional bias (positive at z₀=2), gate-ordering /
        # insert-delete asymmetry is confirmed at L=4 — localization
        # complete. If empirical ⟨N⟩ matches exact, the move kernel is
        # correct even in this regime, and the bug is in NS-loop behavior
        # at the transitional regime not captured by kernel in isolation.
        z0 = 2.0
        ground_U = _l4_ground_state_U(ham, square_lattice)
        u_max = ground_U + 0.1

        exact_mean_N = _exact_constrained_mean_N(z0, u_max, ham, square_lattice)
        exact_PN     = _exact_constrained_N_distribution(z0, u_max, ham, square_lattice)

        n_chains = 50
        n_steps_per_chain = 5000
        Ns = Int[]
        n_sites = num_sites(square_lattice)
        p0 = z0 / (1.0 + z0)

        ins_acc_total = 0
        del_acc_total = 0
        ins_total_estimated = 0.0
        del_total_estimated = 0.0

        Random.seed!(8102)
        for _ in 1:n_chains
            walker = LatticeWalker(deepcopy(square_lattice), energy=0.0u"eV", iter=0)
            random_microstate!(walker.configuration; p=p0)
            assign_energy!(walker, ham)
            attempts = 0
            while walker.energy.val >= u_max && attempts < 50000
                random_microstate!(walker.configuration; p=p0)
                assign_energy!(walker, ham)
                attempts += 1
            end
            @assert walker.energy.val < u_max  "Initial walker exceeds u_max — adjust ceiling"

            for _ in 1:n_steps_per_chain
                N_before = sum(walker.configuration.components[1])
                ins_total_estimated += (n_sites - N_before) / n_sites
                del_total_estimated += N_before / n_sites
                _, _, walker, _, _ = MC_ideal_gas_reference_walk!(
                    1, walker, ham, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                N_after = sum(walker.configuration.components[1])
                if N_after - N_before == 1
                    ins_acc_total += 1
                elseif N_after - N_before == -1
                    del_acc_total += 1
                end
                push!(Ns, N_after)
            end
        end

        emp_PN = zeros(Float64, n_sites + 1)
        for N in Ns
            emp_PN[N + 1] += 1
        end
        emp_PN ./= length(Ns)
        tvd = 0.5 * sum(abs.(emp_PN .- exact_PN))

        ins_acc_rate = ins_total_estimated > 0 ? ins_acc_total / ins_total_estimated : NaN
        del_acc_rate = del_total_estimated > 0 ? del_acc_total / del_total_estimated : NaN

        # Diagnostic logging — read this regardless of pass/fail to inspect
        # the insert-vs-delete asymmetry.
        @info "B2 deep-compression z₀=2" u_max=u_max ground_U=ground_U exact_mean_N=exact_mean_N empirical_mean_N=mean(Ns) signed_diff=mean(Ns) - exact_mean_N tvd=tvd ins_acc_rate=ins_acc_rate del_acc_rate=del_acc_rate ins_acc=ins_acc_total del_acc=del_acc_total ins_est=ins_total_estimated del_est=del_total_estimated

        # Tight equilibrium assertions. With ~250K steps and n_chains=50,
        # statistical noise on ⟨N⟩ is well below 0.05. If gate-ordering
        # introduces a per-iteration bias of order 1e-4 (consistent with
        # the L=6 transitional bias scaling), it would accumulate
        # detectably — failure of these assertions IS the localization
        # signal. Match against analytical equilibrium under the
        # constrained prior with the same tight tolerance as the z₀=1
        # control above.
        @test abs(mean(Ns) - exact_mean_N) ≤ 0.05
        @test tvd ≤ 0.05
    end

    # ================================================================
    @testset "B2 deep-compression: L=6 z₀=1, u_max ≈ ground+0.2 (control)" begin
        # L=6 control for the deep-compression test below. Same u_max as
        # the z₀=2 test, with z₀=1 (Metropolis trivial). Confirms the
        # L=6 kernel-isolated harness is sound: at z₀=1 the kernel should
        # reproduce the analytical conditional regardless of L. If the
        # z₀=2 test below deviates while this control passes, the bias is
        # specific to z₀ ≠ 1 at L=6 deep compression — the prediction of
        # the L-dependent kernel stationary distribution hypothesis.
        L6_lattice = MLattice{1,SquareLattice}(
            lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(6, 6, 1), periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:36]],
            adsorptions=:full)
        n_sites = num_sites(L6_lattice)

        z0 = 1.0
        ground_U = _l6_ground_state_U(ham)
        u_max = ground_U + 0.2

        exact_mean_N, exact_PN = _exact_constrained_distribution_l6(
            z0, u_max, ham, L6_lattice; N_lo=30)

        # Initialize at full lattice (N=36, U=ground); kernel relaxes
        # to conditional within the chain length. Random-microstate
        # rejection-sampling at p=0.5 is too inefficient at L=6 deep
        # compression (acceptance rate ≪ 10⁻⁴ for U<u_max).
        n_chains = 25
        n_chunks = 50
        steps_per_chunk = 100
        burn_in_chunks = 10
        Ns = Int[]

        Random.seed!(8201)
        for _ in 1:n_chains
            walker = LatticeWalker(deepcopy(L6_lattice), energy=0.0u"eV", iter=0)
            walker.configuration.components[1] .= true   # full lattice; U=ground<u_max
            assign_energy!(walker, ham)
            @assert walker.energy.val < u_max

            for chunk in 1:n_chunks
                _, _, walker, _, _ = MC_ideal_gas_reference_walk!(
                    steps_per_chunk, walker, ham, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                if chunk > burn_in_chunks
                    push!(Ns, sum(walker.configuration.components[1]))
                end
            end
        end

        emp_PN = zeros(Float64, n_sites + 1)
        for N in Ns
            emp_PN[N + 1] += 1
        end
        emp_PN ./= length(Ns)
        tvd = 0.5 * sum(abs.(emp_PN .- exact_PN))

        @info "B2 L=6 deep-compression z₀=1 control" ground_U=ground_U u_max=u_max exact_mean_N=exact_mean_N empirical_mean_N=mean(Ns) tvd=tvd n_samples=length(Ns)

        # Generous tolerance — control should pass with margin
        @test abs(mean(Ns) - exact_mean_N) ≤ 0.1
        @test tvd ≤ 0.05
    end

    # ================================================================
    @testset "B2 deep-compression: L=6 z₀=2, u_max ≈ ground+0.2" begin
        # The decisive test for the kernel-asymmetry-under-moving-ceiling
        # vs NS-theory-formula-at-z₀≠1 disambiguation.
        #
        # Setup: u_max = ground+0.2 (≈ −2.14 eV at L=6) — a richer
        # constrained support than L=4's ground+0.1 (which only contained
        # 17 configs at N ∈ {15, 16}). At u_max=−2.14, configs with
        # U < u_max include all N=36 (1 config), N=35 (36 configs),
        # N=34 (630 configs), and a subset of N=33 configs. Multiple N
        # values; non-trivial constrained manifold.
        #
        # If the kernel is correct at L=6 (matches analytical conditional),
        # the bug is in NS-theory's compression-rate formula at z₀≠1
        # (independent of the kernel itself).
        #
        # If the kernel deviates at L=6 (mismatches conditional), the
        # kernel has L-dependent stationary problem — refines the
        # localization to a kernel internal at L=6 deep compression.
        L6_lattice = MLattice{1,SquareLattice}(
            lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(6, 6, 1), periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:36]],
            adsorptions=:full)
        n_sites = num_sites(L6_lattice)

        z0 = 2.0
        ground_U = _l6_ground_state_U(ham)
        u_max = ground_U + 0.2

        exact_mean_N, exact_PN = _exact_constrained_distribution_l6(
            z0, u_max, ham, L6_lattice; N_lo=30)

        n_chains = 25
        n_chunks = 50
        steps_per_chunk = 100
        burn_in_chunks = 10
        Ns = Int[]

        Random.seed!(8202)
        for _ in 1:n_chains
            walker = LatticeWalker(deepcopy(L6_lattice), energy=0.0u"eV", iter=0)
            walker.configuration.components[1] .= true
            assign_energy!(walker, ham)
            @assert walker.energy.val < u_max

            for chunk in 1:n_chunks
                _, _, walker, _, _ = MC_ideal_gas_reference_walk!(
                    steps_per_chunk, walker, ham, u_max, z0;
                    p_flip=1.0, swaps_freq=0, clusters_freq=0,
                    energy_perturb=0.0)
                if chunk > burn_in_chunks
                    push!(Ns, sum(walker.configuration.components[1]))
                end
            end
        end

        emp_PN = zeros(Float64, n_sites + 1)
        for N in Ns
            emp_PN[N + 1] += 1
        end
        emp_PN ./= length(Ns)
        tvd = 0.5 * sum(abs.(emp_PN .- exact_PN))

        @info "B2 L=6 deep-compression z₀=2" ground_U=ground_U u_max=u_max exact_mean_N=exact_mean_N empirical_mean_N=mean(Ns) signed_diff=mean(Ns) - exact_mean_N tvd=tvd n_samples=length(Ns)

        # Same tolerance as the z₀=1 control. If the kernel is L-correct,
        # both pass. If the kernel has L=6 deep-compression bias, this
        # one fails — and the magnitude/direction of failure is the
        # localization signal.
        @test abs(mean(Ns) - exact_mean_N) ≤ 0.1
        @test tvd ≤ 0.05
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
    @testset "_ideal_gas_reference_log_boltzmann_weights helper" begin
        # Hand-computed reference: log w̃_j = log ω_j + N_j(βμ - log z₀) - β U_j.
        β  = 1.0
        μ  = 0.5
        z0 = 2.0
        ωi = [0.5, 0.3, 0.2]
        Us = [0.0, 1.0, -1.0]
        Ns = [0, 1, 2]
        log_z0 = log(z0)
        expected = [log(ωi[i]) + Ns[i] * (β*μ - log_z0) - β*Us[i] for i in 1:3]

        got = AnalysisTools._ideal_gas_reference_log_boltzmann_weights(
            β, μ, ωi, Us, Ns; z0=z0)
        @test got isa Vector{Float64}
        @test length(got) == 3
        @test got ≈ expected rtol=1e-12

        # Validation
        @test_throws DimensionMismatch AnalysisTools._ideal_gas_reference_log_boltzmann_weights(
            β, μ, ωi, [0.0, 1.0], Ns; z0=z0)
        @test_throws ArgumentError AnalysisTools._ideal_gas_reference_log_boltzmann_weights(
            β, μ, ωi, Us, Ns; z0=0.0)
        @test_throws ArgumentError AnalysisTools._ideal_gas_reference_log_boltzmann_weights(
            β, μ, ωi, Us, Ns; z0=-1.0)

        # Empty input → empty result (no error)
        empty_got = AnalysisTools._ideal_gas_reference_log_boltzmann_weights(
            β, μ, Float64[], Float64[], Int[]; z0=z0)
        @test empty_got == Float64[]
    end

    # ================================================================
    @testset "_ideal_gas_reference_log_importance_ratio helper" begin
        # log r_j = N_j(log z - log z₀) - U_j(1/(kbT) - 1/(kbT₀)).
        β   = 1.0 / (8.617333262e-5 * 300.0)  # T = 300 K
        μ   = 0.05
        z0  = 1.5
        T   = 300.0
        T0  = 250.0
        kb  = 8.617333262e-5
        Us  = [-1.0, 0.0, 0.5, 1.0]
        Ns  = [3, 2, 1, 0]
        log_z  = β * μ
        log_z0 = log(z0)
        delta_inv = 1.0/(kb*T) - 1.0/(kb*T0)
        expected = [Ns[i] * (log_z - log_z0) - Us[i] * delta_inv for i in 1:4]

        got = AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, Us, Ns; z0=z0, T=T, T0=T0)
        @test got isa Vector{Float64}
        @test length(got) == 4
        @test got ≈ expected rtol=1e-12

        # T0 default = T → no temperature reweighting term
        got_default = AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, Us, Ns; z0=z0, T=T)
        expected_default = [Ns[i] * (log_z - log_z0) for i in 1:4]
        @test got_default ≈ expected_default rtol=1e-12

        got_explicit = AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, Us, Ns; z0=z0, T=T, T0=T)
        @test got_default ≈ got_explicit rtol=1e-14

        # Validation
        @test_throws DimensionMismatch AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, [0.0, 1.0], Ns; z0=z0, T=T)
        @test_throws ArgumentError AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, Us, Ns; z0=0.0, T=T)
        @test_throws ArgumentError AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, Us, Ns; z0=z0, T=-1.0)
        @test_throws ArgumentError AnalysisTools._ideal_gas_reference_log_importance_ratio(
            β, μ, Us, Ns; z0=z0, T=T, T0=0.0)
    end

    # ================================================================
    @testset "ideal_gas_reference_thermodynamic_stats return_variances scalar" begin
        # The first four scalars must match the default 4-tuple exactly,
        # and the variances must be finite, non-negative, and consistent
        # with var_mean_E ≈ Var_w(U)/N_eff at the sampling distribution.
        β  = 1.0 / (8.617333262e-5 * 300.0)
        μ  = -0.02
        z0 = 1.5
        ωi = [0.4, 0.3, 0.2, 0.1]
        Us = [-1.0, -0.3, 0.2, 0.7]
        Ns = [3, 2, 1, 0]

        out4 = ideal_gas_reference_thermodynamic_stats(
            β, μ, ωi, Us, Ns; z0=z0)
        out7 = ideal_gas_reference_thermodynamic_stats(
            β, μ, ωi, Us, Ns; z0=z0, return_variances=true)

        @test length(out4) == 4
        @test length(out7) == 7
        @test all(isfinite, out7)

        # First four match exactly
        for i in 1:4
            @test out7[i] == out4[i]
        end

        # Variances are non-negative
        @test out7[5] ≥ 0.0   # var_mean_E
        @test out7[6] ≥ 0.0   # var_mean_N
        @test out7[7] ≥ 0.0   # var_ln_Z_relative

        # Empty input → 7 NaNs
        out_empty = ideal_gas_reference_thermodynamic_stats(
            β, μ, Float64[], Float64[], Int[]; z0=z0, return_variances=true)
        @test length(out_empty) == 7
        @test all(isnan, out_empty)

        # Existing 4-tuple behavior preserved when return_variances=false
        out_empty4 = ideal_gas_reference_thermodynamic_stats(
            β, μ, Float64[], Float64[], Int[]; z0=z0)
        @test length(out_empty4) == 4
        @test all(isnan, out_empty4)
    end

    # ================================================================
    @testset "ideal_gas_reference_thermodynamic_stats return_variances DataFrame" begin
        df = DataFrame(
            iter=[1, 2, 3, 4, 5, 6],
            emax=[-1.0, -0.5, 0.0, 0.3, 0.5, 1.0],
            num_particles=[5, 4, 3, 2, 1, 0])

        μs = [-0.02, 0.0, 0.02]
        Ts = [200.0, 300.0]

        result4 = ideal_gas_reference_thermodynamic_stats(
            df, μs, Ts, 10, 1.5)
        result7 = ideal_gas_reference_thermodynamic_stats(
            df, μs, Ts, 10, 1.5; return_variances=true)

        @test result7 isa NamedTuple
        @test haskey(result7, :var_mean_E)
        @test haskey(result7, :var_mean_N)
        @test haskey(result7, :var_ln_Z_relative)
        @test size(result7.var_mean_E) == (3, 2)
        @test size(result7.var_mean_N) == (3, 2)
        @test size(result7.var_ln_Z_relative) == (3, 2)
        @test all(>=(0.0), result7.var_mean_E)
        @test all(>=(0.0), result7.var_mean_N)
        @test all(>=(0.0), result7.var_ln_Z_relative)

        # Default 4-field NamedTuple is unchanged
        @test result4.mean_E == result7.mean_E
        @test result4.mean_N == result7.mean_N
        @test result4.Cv_mu == result7.Cv_mu
        @test result4.ln_Z_relative == result7.ln_Z_relative
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
    @testset "combined_..._stats input validation" begin
        df_ok = DataFrame(iter=[1,2,3], emax=[-1.0, 0.0, 1.0], num_particles=[3,2,1])

        # Length mismatch
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok, df_ok], [1.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:inverse_variance)

        # Empty
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            DataFrame[], Float64[], [-0.05, 0.0], [300.0], 10, 16;
            method=:inverse_variance)

        # Non-positive z0
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [0.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:inverse_variance)
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [-1.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:inverse_variance)

        # Bad n_sites
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [1.0], [-0.05, 0.0], [300.0], 10, 0;
            method=:inverse_variance)
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [1.0], [-0.05, 0.0], [300.0], 10, -3;
            method=:inverse_variance)

        # Bad min_N_eff
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [1.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:inverse_variance, min_N_eff=-0.1)

        # Unknown method
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [1.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:not_a_method)
    end

    # ================================================================
    @testset "Method A: 1-run reduction to single-run estimation" begin
        # With K=1, :inverse_variance must reduce exactly to single-run
        # `ideal_gas_reference_thermodynamic_stats(...; return_variances=true)`.
        df = DataFrame(
            iter=collect(1:30),
            emax=collect(range(-2.0, 2.0, length=30)),
            num_particles=[mod(i, 6) for i in 1:30])

        μs = [-0.04, -0.02, 0.0, 0.02, 0.04]
        Ts = [200.0, 300.0, 500.0]
        n_walkers = 50
        z0 = 1.5

        single = ideal_gas_reference_thermodynamic_stats(
            df, μs, Ts, n_walkers, z0; return_variances=true)
        combined = combined_ideal_gas_reference_thermodynamic_stats(
            [df], [z0], μs, Ts, n_walkers, 16; method=:inverse_variance)

        @test combined.mean_E ≈ single.mean_E rtol=1e-12
        @test combined.mean_N ≈ single.mean_N rtol=1e-12
        @test combined.Cv_mu ≈ single.Cv_mu rtol=1e-12
        @test combined.ln_Z_relative ≈ single.ln_Z_relative rtol=1e-12
        @test combined.var_mean_E ≈ single.var_mean_E rtol=1e-12 atol=1e-30
        @test combined.var_mean_N ≈ single.var_mean_N rtol=1e-12 atol=1e-30
        @test combined.var_ln_Z_relative ≈ single.var_ln_Z_relative rtol=1e-12 atol=1e-30

        # N_eff_per_run shape
        @test size(combined.N_eff_per_run) == (5, 3, 1)
    end

    # ================================================================
    @testset "Method A: 2 identical runs halve the variance" begin
        # Combining a run with itself doubles the precision; combined
        # variance should be (per-run variance)/2 to within fp noise.
        df = DataFrame(
            iter=collect(1:30),
            emax=collect(range(-2.0, 2.0, length=30)),
            num_particles=[mod(i, 6) for i in 1:30])

        μs = [-0.02, 0.0, 0.02]
        Ts = [300.0]
        n_walkers = 50
        z0 = 1.5

        single = ideal_gas_reference_thermodynamic_stats(
            df, μs, Ts, n_walkers, z0; return_variances=true)
        combined = combined_ideal_gas_reference_thermodynamic_stats(
            [df, df], [z0, z0], μs, Ts, n_walkers, 16;
            method=:inverse_variance)

        # Means should equal single-run means (averaging identical runs).
        @test combined.mean_E ≈ single.mean_E rtol=1e-12
        @test combined.mean_N ≈ single.mean_N rtol=1e-12
        @test combined.ln_Z_relative ≈ single.ln_Z_relative rtol=1e-12

        # Variances should be halved (or zero if single is zero).
        for i in eachindex(single.var_mean_E)
            if single.var_mean_E[i] > 0
                @test combined.var_mean_E[i] ≈ single.var_mean_E[i] / 2 rtol=1e-12
            else
                @test combined.var_mean_E[i] == 0.0
            end
        end
        for i in eachindex(single.var_mean_N)
            if single.var_mean_N[i] > 0
                @test combined.var_mean_N[i] ≈ single.var_mean_N[i] / 2 rtol=1e-12
            else
                @test combined.var_mean_N[i] == 0.0
            end
        end
    end

    # ================================================================
    @testset "Method A: 2-run combination smoke" begin
        # Two distinct runs at different z₀: combined output is finite
        # and combined variance is at most the smallest per-run variance.
        df1 = DataFrame(
            iter=collect(1:30),
            emax=collect(range(-2.0, 2.0, length=30)),
            num_particles=[mod(i, 6) for i in 1:30])
        df2 = DataFrame(
            iter=collect(1:30),
            emax=collect(range(-1.5, 2.5, length=30)),
            num_particles=[mod(i + 1, 6) for i in 1:30])

        μs = [-0.02, 0.0, 0.02]
        Ts = [300.0]
        n_walkers = 50
        z0s = [1.0, 2.0]

        s1 = ideal_gas_reference_thermodynamic_stats(
            df1, μs, Ts, n_walkers, z0s[1]; return_variances=true)
        s2 = ideal_gas_reference_thermodynamic_stats(
            df2, μs, Ts, n_walkers, z0s[2]; return_variances=true)

        combined = combined_ideal_gas_reference_thermodynamic_stats(
            [df1, df2], z0s, μs, Ts, n_walkers, 16;
            method=:inverse_variance)

        # Output structure
        @test combined isa NamedTuple
        @test haskey(combined, :mean_E)
        @test haskey(combined, :var_mean_E)
        @test haskey(combined, :N_eff_per_run)
        @test size(combined.mean_E) == (3, 1)
        @test size(combined.N_eff_per_run) == (3, 1, 2)

        # All finite
        @test all(isfinite, combined.mean_E)
        @test all(isfinite, combined.mean_N)
        @test all(isfinite, combined.ln_Z_relative)

        # Combined variance at most the smallest per-run variance
        for i in eachindex(combined.var_mean_E)
            min_var = min(s1.var_mean_E[i], s2.var_mean_E[i])
            if min_var > 0  # well-defined comparison
                @test combined.var_mean_E[i] ≤ min_var + 1e-14
            end
        end
        for i in eachindex(combined.var_mean_N)
            min_var = min(s1.var_mean_N[i], s2.var_mean_N[i])
            if min_var > 0
                @test combined.var_mean_N[i] ≤ min_var + 1e-14
            end
        end
    end

    # ================================================================
    @testset "Method A: return_per_run=true exposes per-run estimates" begin
        df1 = DataFrame(
            iter=collect(1:20), emax=collect(range(-1.0, 1.0, length=20)),
            num_particles=[mod(i, 5) for i in 1:20])
        df2 = DataFrame(
            iter=collect(1:20), emax=collect(range(-0.5, 1.5, length=20)),
            num_particles=[mod(i, 5) for i in 1:20])

        result = combined_ideal_gas_reference_thermodynamic_stats(
            [df1, df2], [1.0, 2.0], [-0.02, 0.02], [300.0], 50, 16;
            method=:inverse_variance, return_per_run=true)

        @test haskey(result, :per_run)
        @test length(result.per_run) == 2
        for k in 1:2
            @test result.per_run[k] isa NamedTuple
            @test haskey(result.per_run[k], :mean_E)
            @test haskey(result.per_run[k], :var_mean_E)
        end

        # Without return_per_run=true, the field is absent
        result_no_pr = combined_ideal_gas_reference_thermodynamic_stats(
            [df1, df2], [1.0, 2.0], [-0.02, 0.02], [300.0], 50, 16;
            method=:inverse_variance, return_per_run=false)
        @test !haskey(result_no_pr, :per_run)
    end

    # ================================================================
    @testset "Method A: min_N_eff filtering" begin
        # When min_N_eff is set very high, no run survives at any cell;
        # all results should be NaN.
        df = DataFrame(
            iter=collect(1:20), emax=collect(range(-1.0, 1.0, length=20)),
            num_particles=[mod(i, 5) for i in 1:20])

        # No filtering — finite outputs
        result_no = combined_ideal_gas_reference_thermodynamic_stats(
            [df, df], [1.0, 2.0], [-0.02, 0.02], [300.0], 50, 16;
            method=:inverse_variance, min_N_eff=0.0)
        @test all(isfinite, result_no.mean_E)

        # Aggressive filter — every run rejected, every cell NaN
        result_strict = combined_ideal_gas_reference_thermodynamic_stats(
            [df, df], [1.0, 2.0], [-0.02, 0.02], [300.0], 50, 16;
            method=:inverse_variance, min_N_eff=1e10)
        @test all(isnan, result_strict.mean_E)
        @test all(isnan, result_strict.mean_N)
        @test all(isnan, result_strict.ln_Z_relative)
        @test all(isnan, result_strict.var_mean_E)
    end

    # ================================================================
    @testset "Method B: K=1 reduces to single-run absolute Ξ" begin
        # With K=1 and α_1=N_total, the mixture reduces to π_0^(1):
        #   log f_mix^unn(σ) = N(σ) · log z₀^(1)
        #   log Z_mix        = M · log(1 + z₀^(1))
        # so the combined absolute Ξ must equal the single-run formula
        # M·log(1+z₀) + ln_Z_relative.
        df = DataFrame(
            iter=collect(1:30),
            emax=collect(range(-2.0, 2.0, length=30)),
            num_particles=[mod(i, 6) for i in 1:30])
        z0 = 1.5
        M = 16
        n_walkers = 50
        μs = [-0.04, -0.02, 0.0, 0.02, 0.04]
        Ts = [200.0, 300.0, 500.0]

        single = ideal_gas_reference_thermodynamic_stats(
            df, μs, Ts, n_walkers, z0)
        combined = combined_ideal_gas_reference_thermodynamic_stats(
            [df], [z0], μs, Ts, n_walkers, M; method=:mixture_importance)

        single_ln_Ξ = M * log(1.0 + z0) .+ single.ln_Z_relative

        @test combined.mean_E ≈ single.mean_E rtol=1e-12
        @test combined.mean_N ≈ single.mean_N rtol=1e-12
        @test combined.Cv_mu ≈ single.Cv_mu rtol=1e-10 atol=1e-15
        @test combined.ln_Z_relative ≈ single.ln_Z_relative rtol=1e-12
        @test combined.ln_Z_absolute ≈ single_ln_Ξ rtol=1e-12

        # N_eff_pooled at K=1 ≈ N_eff^(1) (same per-sample weights)
        @test size(combined.N_eff_pooled) == (5, 3)
        @test all(isfinite, combined.N_eff_pooled)
        @test all(>(0), combined.N_eff_pooled)
    end

    # ================================================================
    @testset "Method B: K=1 reduction at L=8 scale (regression guard)" begin
        # The original K=1 reduction test above uses M=16, n_walkers=50,
        # 30 samples — small enough that any algorithmic divergence between
        # the two paths would show up at <1 ulp. This guard exercises the
        # same identity at L=8 scale (M=64, n_walkers=100, ~1685 samples)
        # over the (μ, T) grid used by the head-to-head benchmark in
        # `scripts/benchmark_gcns_combined_vs_omega.jl`, where β·U dynamic
        # range reaches ~150 and ln(ω) spans ~20 orders of magnitude.
        # Synthetic but deterministic data; no NS run required.
        M = 64
        n_walkers = 100
        n_samples = 1685
        iter = collect(1:n_samples)
        # Deterministic N(i) sweeping the full [0, M] range many times.
        num_particles = [Int(mod(7*i + 23 + i*(i-1) ÷ 5, M + 1))
                         for i in 1:n_samples]
        # Realistic-magnitude energies: roughly −ε_ads·N − pair·N²/M plus a
        # small deterministic ripple. |U|max ≈ 2.7 eV; at T=200 K this gives
        # |β·U| ≈ 156, the regime that motivated this guard.
        emax = [-0.04 * num_particles[i] - 0.0025 * num_particles[i]^2 / M +
                0.003 * (mod(i * 13 + 7, 101) / 101.0 - 0.5)
                for i in 1:n_samples]
        df = DataFrame(iter=iter, emax=emax, num_particles=num_particles)

        μs = [-0.06, -0.04, -0.02, 0.0, 0.02]
        Ts = [200.0, 300.0, 400.0]

        # Sweep z₀ at and off-center: at z₀=1 the fugacity-ratio term
        # vanishes, so off-center cases are the ones that exercise
        # `N·(β·μ − log z₀)` against `−log π_mix^norm` with a non-trivial
        # `log z₀` contribution.
        for z0 in (1.0, 0.5, 2.0)
            single = ideal_gas_reference_thermodynamic_stats(
                df, μs, Ts, n_walkers, z0)
            combined = combined_ideal_gas_reference_thermodynamic_stats(
                [df], [z0], μs, Ts, n_walkers, M; method=:mixture_importance)

            @test combined.mean_E ≈ single.mean_E rtol=1e-12
            @test combined.mean_N ≈ single.mean_N rtol=1e-12
            @test combined.Cv_mu ≈ single.Cv_mu rtol=1e-10 atol=1e-15
            @test combined.ln_Z_relative ≈ single.ln_Z_relative rtol=1e-12
            # ln_Z bookkeeping: single's relative + M·log(1+z₀) is the
            # combined path's absolute Ξ at K=1. Documented in the
            # `_combine_mixture_importance` docstring.
            @test combined.ln_Z_absolute ≈
                  single.ln_Z_relative .+ M * log1p(z0) rtol=1e-12
        end
    end

    # ================================================================
    @testset "Method B: K=2 same df, same z₀ reduces to K=1" begin
        # Combining a run with itself at the same z₀ gives the same
        # absolute Ξ as a single-run estimate at that z₀:
        #   log_alpha = log(1/2) cancels against the same factor in
        #   log Z_mix = log(2) + log(1/2) + M·log(1+z₀) = M·log(1+z₀).
        # Per-sample log w̃ contains an extra log(1/2) but each sample
        # appears in both copies of the df, restoring the K=1 sum.
        df = DataFrame(
            iter=collect(1:30),
            emax=collect(range(-2.0, 2.0, length=30)),
            num_particles=[mod(i, 6) for i in 1:30])
        z0 = 2.0
        M = 16
        n_walkers = 50
        μs = [-0.02, 0.0, 0.02]
        Ts = [300.0]

        single = combined_ideal_gas_reference_thermodynamic_stats(
            [df], [z0], μs, Ts, n_walkers, M; method=:mixture_importance)
        doubled = combined_ideal_gas_reference_thermodynamic_stats(
            [df, df], [z0, z0], μs, Ts, n_walkers, M; method=:mixture_importance)

        @test doubled.mean_E ≈ single.mean_E rtol=1e-12
        @test doubled.mean_N ≈ single.mean_N rtol=1e-12
        @test doubled.Cv_mu ≈ single.Cv_mu rtol=1e-10 atol=1e-15
        @test doubled.ln_Z_relative ≈ single.ln_Z_relative rtol=1e-12
        @test doubled.ln_Z_absolute ≈ single.ln_Z_absolute rtol=1e-12
    end

    # ================================================================
    @testset "Method B: K=2 different z₀ — analytic Z_mix constant" begin
        # Cross-check the absolute-Ξ additive constant against direct
        # evaluation of LSE_k[log(α_k/N_total) + M·log(1+z₀^(k))].
        # We can't easily check the sample-level estimator without running
        # NS, but we can verify the additive constant is wired correctly
        # by examining ln_Z_absolute - ln_Z_relative.
        df1 = DataFrame(
            iter=collect(1:20), emax=collect(range(-1.0, 1.0, length=20)),
            num_particles=[mod(i, 5) for i in 1:20])
        df2 = DataFrame(
            iter=collect(1:30), emax=collect(range(-1.5, 1.5, length=30)),
            num_particles=[mod(i, 5) for i in 1:30])

        z0s = [0.5, 2.0]
        M = 16
        n_walkers = 50
        μs = [0.0]
        Ts = [300.0]

        result = combined_ideal_gas_reference_thermodynamic_stats(
            [df1, df2], z0s, μs, Ts, n_walkers, M; method=:mixture_importance)

        # Hand-compute the additive constant
        αs = Float64[20, 30]
        N_total = sum(αs)
        log_alpha_over_total = log.(αs ./ N_total)
        log_Z_mix_terms = [log_alpha_over_total[k] + M * log(1.0 + z0s[k]) for k in 1:2]
        m = maximum(log_Z_mix_terms)
        log_Z_mix = m + log(sum(exp.(log_Z_mix_terms .- m)))

        @test result.ln_Z_absolute[1, 1] ≈ result.ln_Z_relative[1, 1] + log_Z_mix rtol=1e-12
    end

    # ================================================================
    @testset "Method B: input validation parity with Method A" begin
        df_ok = DataFrame(iter=[1,2,3], emax=[-1.0, 0.0, 1.0], num_particles=[3,2,1])
        # length mismatch
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok, df_ok], [1.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:mixture_importance)
        # bad z0
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [-1.0], [-0.05, 0.0], [300.0], 10, 16;
            method=:mixture_importance)
        # bad n_sites
        @test_throws ArgumentError combined_ideal_gas_reference_thermodynamic_stats(
            [df_ok], [1.0], [-0.05, 0.0], [300.0], 10, 0;
            method=:mixture_importance)
    end

    # ================================================================
    @testset "Method B: return_per_run=true" begin
        df1 = DataFrame(
            iter=collect(1:20), emax=collect(range(-1.0, 1.0, length=20)),
            num_particles=[mod(i, 5) for i in 1:20])
        df2 = DataFrame(
            iter=collect(1:20), emax=collect(range(-0.5, 1.5, length=20)),
            num_particles=[mod(i, 5) for i in 1:20])

        result = combined_ideal_gas_reference_thermodynamic_stats(
            [df1, df2], [1.0, 2.0], [-0.02, 0.02], [300.0], 50, 16;
            method=:mixture_importance, return_per_run=true)

        @test haskey(result, :ln_Z_absolute)
        @test haskey(result, :N_eff_pooled)
        @test haskey(result, :per_run)
        @test length(result.per_run) == 2

        result_no_pr = combined_ideal_gas_reference_thermodynamic_stats(
            [df1, df2], [1.0, 2.0], [-0.02, 0.02], [300.0], 50, 16;
            method=:mixture_importance, return_per_run=false)
        @test !haskey(result_no_pr, :per_run)
    end

    # ================================================================
    @testset "_logsumexp helper" begin
        lse = AnalysisTools._logsumexp

        # Basic correctness against direct evaluation
        @test lse([0.0]) == 0.0
        @test lse([0.0, 0.0]) ≈ log(2.0)
        @test lse([0.0, 0.0, 0.0]) ≈ log(3.0)
        @test lse([log(1.0), log(2.0), log(3.0)]) ≈ log(6.0)
        @test lse([log(0.1), log(0.5), log(2.0)]) ≈ log(0.1 + 0.5 + 2.0)

        # Numerical stability under large shifts
        @test lse([1000.0, 1000.0]) ≈ 1000.0 + log(2.0)
        @test lse([-1000.0, -1000.0]) ≈ -1000.0 + log(2.0)
        # Direct exp would overflow; LSE must not.
        @test isfinite(lse([800.0, 700.0, 750.0]))
        @test lse([800.0, 700.0, 750.0]) ≈ 800.0 + log(1.0 + exp(-100.0) + exp(-50.0))

        # Heterogeneous magnitudes
        @test lse([1e3, -1e3]) ≈ 1e3   # the large one dominates
        @test lse([1e3, -1e3]) - 1e3 < 1e-300

        # Edge cases
        @test lse(Float64[]) == -Inf
        @test lse([-Inf, -Inf]) == -Inf
        @test lse([-Inf, 0.0]) == 0.0   # one finite element
        @test lse([+Inf, 0.0]) == +Inf

        # NaN propagation
        @test isnan(lse([0.0, NaN]))
    end

    # ================================================================
    @testset "_log_mixture_prior_density helper" begin
        lmp = AnalysisTools._log_mixture_prior_density

        # K=1: log f_mix^unn(σ_j) = log(α_1/N_total) + N_j · log(z₀^(1))
        # With α_1 = N_total, log(α_1/N_total) = 0.
        log_alpha = [0.0]
        log_z0s   = [log(1.5)]
        for N_j in 0:6
            @test lmp(N_j, log_alpha, log_z0s) ≈ N_j * log(1.5) atol=1e-14
        end

        # K=2 with hand-computed reference
        log_alpha2 = [log(0.4), log(0.6)]
        log_z0s2   = [log(1.0), log(2.0)]
        for N_j in 0:5
            expected = log(0.4 * 1.0^N_j + 0.6 * 2.0^N_j)
            @test lmp(N_j, log_alpha2, log_z0s2) ≈ expected atol=1e-14
        end

        # K=3 asymmetric: weights 0.5, 0.3, 0.2 at z₀ = 0.5, 1.0, 4.0
        log_alpha3 = [log(0.5), log(0.3), log(0.2)]
        log_z0s3   = [log(0.5), log(1.0), log(4.0)]
        for N_j in 0:8
            expected = log(0.5 * 0.5^N_j + 0.3 * 1.0^N_j + 0.2 * 4.0^N_j)
            @test lmp(N_j, log_alpha3, log_z0s3) ≈ expected atol=1e-14
        end

        # Consistency: equal αs, equal z₀s ⇒ log(K · 1/K · z₀^N) = N·log(z₀)
        K = 4
        log_alpha_eq = fill(-log(K), K)   # log(α_k/N_total) = log(1/K)
        log_z0s_eq = fill(log(1.5), K)
        for N_j in 0:8
            @test lmp(N_j, log_alpha_eq, log_z0s_eq) ≈ N_j * log(1.5) atol=1e-14
        end

        # Numerical stability with large N_j and large fugacity span
        log_alpha_big = [log(0.5), log(0.5)]
        log_z0s_big   = [log(1e-2), log(1e2)]
        N_j = 30
        # exp(N_j · log(1e2)) = 10^60 — would overflow naive exp.
        got = lmp(N_j, log_alpha_big, log_z0s_big)
        @test isfinite(got)
        # The N_j=30 term at z₀=100 dominates; result ≈ log(0.5) + 30·log(100).
        @test got ≈ log(0.5) + 30 * log(1e2) atol=1e-10
    end

    # ================================================================
    @testset "_log_mixture_prior_density: K=2 distinct z₀, normalized form" begin
        # The K≥2 analog of the K=1 reduction check. With distinct z₀,
        # the mixture must be evaluated in NORMALIZED form
        #   log π_mix^norm(σ) = LSE_k[log(α_k/N_total) − M·log(1+z₀^(k))
        #                              + N(σ)·log z₀^(k)]
        # to give the unbiased Method B importance weight. This test
        # passes the *normalized* weights to the helper (which then
        # produces log π_mix^norm) and compares against direct
        # hand-computed values at small M and small K=2.
        #
        # An earlier campaign bug used the *unnormalized* form
        # (omitting −M·log(1+z₀)) and produced incorrect ln Ξ for
        # K≥2 with distinct z₀. This test would have caught it.
        lmp = AnalysisTools._log_mixture_prior_density

        M = 4
        α = [1.0, 1.0]            # equal sample counts
        N_total = sum(α)
        z0s = [1.0, 3.0]
        log_alpha_over_total = log.(α ./ N_total)             # both = log(1/2)
        log_normed_weights = [log_alpha_over_total[k] - M * log1p(z0s[k])
                              for k in 1:2]
        log_z0s = log.(z0s)

        # Hand-computed:
        #   π_mix^norm(σ) = (1/2)·1^N/(1+1)^4 + (1/2)·3^N/(1+3)^4
        #                 = (1/2)/16 + (1/2)·3^N/256 = (16 + 3^N)/512
        for N_j in 0:6
            expected = log((16.0 + 3.0^N_j) / 512.0)
            got = lmp(N_j, log_normed_weights, log_z0s)
            @test got ≈ expected atol=1e-12
        end

        # Asymmetric α, well-separated z₀: α=[100, 50], z₀=[0.25, 4.0], M=3.
        # log π_mix^norm(σ) at given N(σ) should equal hand-computed form.
        M2 = 3
        α2 = [100.0, 50.0]
        N_total2 = sum(α2)
        z0s2 = [0.25, 4.0]
        log_alpha_over_total2 = log.(α2 ./ N_total2)        # [log(2/3), log(1/3)]
        log_normed_weights2 = [log_alpha_over_total2[k] - M2 * log1p(z0s2[k])
                               for k in 1:2]
        log_z0s2 = log.(z0s2)

        # π_mix^norm(σ) = (2/3)·(0.25)^N/(1.25)^3 + (1/3)·(4)^N/(5)^3
        for N_j in 0:5
            expected = log((2.0/3.0) * 0.25^N_j / 1.25^3
                          + (1.0/3.0) * 4.0^N_j / 5.0^3)
            got = lmp(N_j, log_normed_weights2, log_z0s2)
            @test got ≈ expected atol=1e-12
        end

        # K=1 reduction of the normalized form: log π_mix^norm at K=1 with
        # α_1 = N_total reduces to N·log z₀ − M·log(1+z₀), the per-site
        # Bernoulli prior log-density.
        for z0 in (0.5, 1.0, 2.5)
            log_normed_weights_K1 = [0.0 - M * log1p(z0)]
            log_z0s_K1 = [log(z0)]
            for N_j in 0:M
                got = lmp(N_j, log_normed_weights_K1, log_z0s_K1)
                expected = N_j * log(z0) - M * log1p(z0)
                @test got ≈ expected atol=1e-14
            end
        end
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
    @testset "Tier 0 multi-run: K=3 z₀ sweep, absolute Ξ formula" begin
        # The most sensitive test of Method B's mixture-aware additive
        # constant. Three NS runs at z₀ ∈ {0.3, 1.0, 3.0} on the
        # non-interacting Hamiltonian; combined absolute Ξ must match
        # the closed-form (1+z(μ,T))^M over a μ window wider than any
        # single run can trust.
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

        z0s = [0.3, 1.0, 3.0]
        n_walkers = 100
        n_ns_steps = Int64(1500)
        T = 300.0
        kb = 8.617333262e-5
        β = 1.0 / (kb * T)

        # Run NS at each z₀
        dfs = DataFrame[]
        for z0 in z0s
            walkers = [LatticeWalker(deepcopy(lattice_template),
                                     energy=0.0u"eV", iter=0)
                       for _ in 1:n_walkers]
            liveset = LatticeGasWalkers(walkers, ham_zero; assign_energy=false)
            params = LatticeGCNSIdealGasReferenceParameters(
                mc_steps=80, z0=z0, energy_perturbation=1e-12,
                allowed_fail_count=10000)
            mc_routine = MCIdealGasReferenceMoves(
                p_flip=1.0, swaps_freq=0, clusters_freq=0)
            save_strat = SaveEveryN(
                df_filename="test_t0mr_$(z0)_df.csv",
                wk_filename="test_t0mr_$(z0).traj",
                ls_filename="test_t0mr_$(z0).ls",
                n_traj=10000, n_snap=10000, n_info=10000)

            df, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
                liveset, params, n_ns_steps, mc_routine, save_strat)
            push!(dfs, df)
        end

        # Wider μ window than any single run can trust: each run's
        # high-N_eff window is roughly μ₀^(k) ± ~0.025 eV, so the
        # combined window of [-0.06, +0.06] requires multiple runs.
        μs = collect(range(-0.06, 0.06, length=7))
        result_B = combined_ideal_gas_reference_thermodynamic_stats(
            dfs, z0s, μs, [T], n_walkers, n_sites;
            method=:mixture_importance)

        # ⟨E⟩ = 0 exactly (no interactions)
        @test all(abs.(result_B.mean_E) .< 1e-6)
        # Cv_μ = 0 exactly (Var(U) = 0)
        @test all(abs.(result_B.Cv_mu) .< 1e-6)

        # ⟨N⟩ matches M·z/(1+z). Tolerance: combined run has N_eff
        # of order K · K_walkers, so SD ≤ ~0.5 at the extremes.
        for (i, μ) in enumerate(μs)
            z = exp(β * μ)
            expected_N = n_sites * z / (1.0 + z)
            @test result_B.mean_N[i, 1] ≈ expected_N atol=1.5
        end

        # *** Critical test: ln_Z_absolute matches M·log(1+z(μ,T)). ***
        # This validates the mixture-aware additive constant
        # log Z_mix = LSE_k[log(α_k/N_total) + M·log(1+z₀^(k))].
        # For non-interacting systems this is exact, so any disagreement
        # beyond NS statistical noise indicates a bug in the additive
        # constant derivation.
        for (i, μ) in enumerate(μs)
            z = exp(β * μ)
            expected_lnΞ = n_sites * log(1.0 + z)
            @test result_B.ln_Z_absolute[i, 1] ≈ expected_lnΞ atol=0.5
        end

        # Method A also produces sensible mean_E, mean_N (no
        # ln_Z_absolute). Method A is fundamentally limited at the
        # μ-window edges: it inverse-variance combines per-run estimates
        # that are biased outside the run's own high-N_eff window
        # (≈ μ₀^(k) ± 0.025 for k = 1..K), and inverse-variance weighting
        # of biased estimators does not recover the truth — a low-variance
        # but biased run can dominate. The 3 z₀ here have μ₀ at
        # ≈ {−0.031, 0, +0.031}, and these windows do not overlap
        # simultaneously, so there is no μ where all 3 runs are
        # well-targeted at once. Method B's mixture-importance estimator
        # resolves this (see test above).
        # Mean_E remains zero exactly under Method A (no interactions).
        # For mean_N we restrict the comparison to the narrow central
        # window [−0.02, +0.02] where all three single-run estimates are
        # reasonably reliable.
        result_A = combined_ideal_gas_reference_thermodynamic_stats(
            dfs, z0s, μs, [T], n_walkers, n_sites;
            method=:inverse_variance)
        @test all(abs.(result_A.mean_E) .< 1e-6)
        method_A_window = findall(μ -> -0.02 - 1e-12 ≤ μ ≤ 0.02 + 1e-12, μs)
        for i in method_A_window
            μ = μs[i]
            z = exp(β * μ)
            expected_N = n_sites * z / (1.0 + z)
            @test result_A.mean_N[i, 1] ≈ expected_N atol=2.0
        end

        # Pooled N_eff should be positive and finite at every cell.
        # The stronger claim — pooled N_eff fills gaps between per-run
        # windows and tends to be ≥ the largest per-run N_eff — is a
        # demonstration target for the Tier 3 benchmark plot, not a
        # strict mathematical guarantee in finite samples (the
        # Veach-Guibas variance-minimization result applies to moment
        # estimators, not directly to ESS).
        @test all(isfinite, result_B.N_eff_pooled)
        @test all(>(0), result_B.N_eff_pooled)

        for z0 in z0s
            rm("test_t0mr_$(z0)_df.csv", force=true)
            rm("test_t0mr_$(z0).traj", force=true)
            rm("test_t0mr_$(z0).ls", force=true)
        end
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
    @testset "Tier 1 multi-run: K=3 z₀ sweep vs exact enumeration" begin
        # 4×4 interacting lattice gas (standard FreeBird GCNS validation
        # Hamiltonian) with three NS runs at z₀ ∈ {0.5, 1.0, 2.0}.
        # Combined estimates from both Method A (inverse variance) and
        # Method B (mixture importance) are compared against exact
        # enumeration over 2^16 = 65536 microstates over a μ-window
        # wider than any single run's high-N_eff region.
        L = 4
        n_sites = L * L
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        lattice_template = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:n_sites]],
            adsorptions=:full)

        T  = 300.0
        kb = 8.617333262e-5
        β  = 1.0 / (kb * T)
        z0s = [0.5, 1.0, 2.0]
        n_walkers = 80
        n_ns_steps = Int64(2000)

        # ----- Exact enumeration -----
        E_states = Vector{Float64}(undef, 2^n_sites)
        N_states = Vector{Int}(undef, 2^n_sites)
        for mask in 0:(2^n_sites - 1)
            lat = deepcopy(lattice_template)
            for site in 1:n_sites
                lat.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_states[mask + 1] = interacting_energy(lat, ham).val
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

        # ----- NS runs at each z₀ -----
        dfs = DataFrame[]
        for z0 in z0s
            walkers = [LatticeWalker(deepcopy(lattice_template),
                                     energy=0.0u"eV", iter=0)
                       for _ in 1:n_walkers]
            liveset = LatticeGasWalkers(walkers, ham; assign_energy=false)
            params = LatticeGCNSIdealGasReferenceParameters(
                mc_steps=80, z0=z0, energy_perturbation=1e-12,
                allowed_fail_count=10000)
            mc_routine = MCIdealGasReferenceMoves(
                p_flip=0.5, swaps_freq=1, clusters_freq=0)
            save_strat = SaveEveryN(
                df_filename="test_t1mr_$(z0)_df.csv",
                wk_filename="test_t1mr_$(z0).traj",
                ls_filename="test_t1mr_$(z0).ls",
                n_traj=10000, n_snap=10000, n_info=10000)
            df, _, _ = ideal_gas_reference_grand_canonical_nested_sampling(
                liveset, params, n_ns_steps, mc_routine, save_strat)
            push!(dfs, df)
        end

        # μ-window wider than any single run's high-N_eff region
        # (single run at z₀=1 has its window roughly μ₀±0.025 = [−0.025, +0.025];
        # combined window with K=3 ≈ [−0.05, +0.05]).
        μs = [-0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06]

        result_B = combined_ideal_gas_reference_thermodynamic_stats(
            dfs, z0s, μs, [T], n_walkers, n_sites;
            method=:mixture_importance)
        result_A = combined_ideal_gas_reference_thermodynamic_stats(
            dfs, z0s, μs, [T], n_walkers, n_sites;
            method=:inverse_variance)

        # ----- Method B vs exact, full wide window -----
        # Headline result: Method B tracks exact enumeration over the wide
        # combined μ-window. Tolerances mirror the single-run Tier 1 atol
        # for ⟨E⟩, ⟨N⟩, and ln Ξ at the center, with slight loosening at
        # the wide-window edges.
        for (i, μ) in enumerate(μs)
            ex_E, ex_N, ex_Cv, ex_lnΞ = exact_stats(μ)

            @test result_B.mean_N[i, 1] ≈ ex_N rtol=0.3 atol=1.5
            @test result_B.mean_E[i, 1] ≈ ex_E rtol=0.3 atol=0.1

            # Cv hardest to converge — sign-check only.
            if isfinite(result_B.Cv_mu[i, 1]) && isfinite(ex_Cv) && ex_Cv > 0
                @test result_B.Cv_mu[i, 1] > 0
            end

            # ln Ξ_absolute is the strongest test of the mixture-aware
            # construction.
            @test result_B.ln_Z_absolute[i, 1] ≈ ex_lnΞ atol=1.0
        end

        # ----- Method A vs exact, narrow central window -----
        # Method A inverse-variance combines per-run estimates and is
        # biased outside individual runs' high-N_eff windows (see Tier 0
        # multi-run docstring). Restrict the comparison to where at
        # least the central run (z₀=1, μ₀=0) provides reliable per-run
        # estimates: μ ∈ [−0.025, +0.025], i.e., |μ| ≤ 0.025.
        method_A_window = findall(μ -> -0.025 - 1e-12 ≤ μ ≤ 0.025 + 1e-12, μs)
        for i in method_A_window
            μ = μs[i]
            ex_E, ex_N, ex_Cv, ex_lnΞ = exact_stats(μ)
            @test result_A.mean_N[i, 1] ≈ ex_N rtol=0.3 atol=2.0
            @test result_A.mean_E[i, 1] ≈ ex_E rtol=0.3 atol=0.1
        end

        # ----- Coverage: combined estimate must be at least as accurate
        # as the single best single-run estimate at every μ.
        # Compare Method B's |error| against the *minimum* per-run |error|
        # (i.e., the closest single run wins each cell). Method B should
        # be at least competitive everywhere.
        per_run_results = [
            ideal_gas_reference_thermodynamic_stats(dfs[k], μs, [T], n_walkers, z0s[k])
            for k in 1:length(z0s)]
        for (i, μ) in enumerate(μs)
            _, ex_N, _, _ = exact_stats(μ)
            err_B = abs(result_B.mean_N[i, 1] - ex_N)
            err_per_run = [abs(per_run_results[k].mean_N[i, 1] - ex_N)
                           for k in 1:length(z0s)]
            best_per_run = minimum(err_per_run)
            # Method B is allowed up to 1.5 absolute units worse than the
            # best per-run estimate (statistical noise margin); in practice
            # it's typically better than the best single run on average.
            @test err_B ≤ best_per_run + 1.5
        end

        for z0 in z0s
            rm("test_t1mr_$(z0)_df.csv", force=true)
            rm("test_t1mr_$(z0).traj", force=true)
            rm("test_t1mr_$(z0).ls", force=true)
        end
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
