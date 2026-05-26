@testset "Atomistic GC-NS, fixed-N construction (ideal gas)" begin

    @testset "IdealGasParameters dispatch" begin
        pot = IdealGasParameters()

        # pair_energy is identically zero
        @test pair_energy(1.0u"Å", pot) == 0.0u"eV"
        @test pair_energy(0.5u"Å", pot) == 0.0u"eV"
        @test pair_energy(100.0u"Å", pot) == 0.0u"eV"

        # Construct a 3-atom periodic system and verify every pairwise dispatch
        # path returns 0.0u"eV".
        box = [[10.0u"Å", 0u"Å", 0u"Å"],
               [0u"Å", 10.0u"Å", 0u"Å"],
               [0u"Å", 0u"Å", 10.0u"Å"]]
        coor = [:Ar => [0.1, 0.2, 0.3],
                :Ar => [0.4, 0.5, 0.6],
                :Ar => [0.7, 0.8, 0.9]]
        sys = FastSystem(periodic_system(coor, box, fractional=true))

        @test interacting_energy(sys, pot) == 0.0u"eV"
        @test interacting_energy(sys, pot, [3], [false]) == 0.0u"eV"
        @test frozen_energy(sys, pot, [3], [true]) == 0.0u"eV"
    end

    @testset "_thermal_wavelength" begin
        # Reference value for H (m=1u) at T=300K: Λ ≈ 1.008 Å.
        Λ_H = FreeBird.AnalysisTools._thermal_wavelength(1.0u"u", 300.0u"K")
        @test isapprox(ustrip(u"Å", Λ_H), 1.008, rtol=1e-2)

        # Exact scaling: Λ ∝ 1/sqrt(m), Λ ∝ 1/sqrt(T).
        Λ_2m = FreeBird.AnalysisTools._thermal_wavelength(2.0u"u", 300.0u"K")
        Λ_2T = FreeBird.AnalysisTools._thermal_wavelength(1.0u"u", 600.0u"K")
        @test isapprox(ustrip(u"Å", Λ_2m) / ustrip(u"Å", Λ_H), 1 / sqrt(2), rtol=1e-10)
        @test isapprox(ustrip(u"Å", Λ_2T) / ustrip(u"Å", Λ_H), 1 / sqrt(2), rtol=1e-10)
    end

    @testset "ideal-gas closed form: Ξ, ⟨N⟩, Var(N), ⟨U⟩" begin
        # Synthesize per-N canonical NS DataFrames for an ideal gas (E ≡ 0).
        # With ω0 = (K+1)/K and many iterations, sum(ω_i) → 1 exactly,
        # making Z_NS^{(N)} = 1 for every N (the analytical answer).
        # N_max = 20 keeps the truncation tail (zV)^N/N! negligible for
        # ⟨N⟩ ≤ 3 at rtol = 1e-3.
        N_max = 20
        N_values = collect(0:N_max)
        K = 120
        n_iters = 5000
        ω0_test = (K + 1) / K

        ns_outputs = [DataFrame(iter=collect(1:n_iters), emax=zeros(n_iters))
                      for _ in N_values]

        V = 1000.0u"Å^3"
        m = 40.0u"u"
        T = 300.0u"K"
        T_grid = [T]

        kb = 8.617333262e-5
        β = 1.0 / (kb * ustrip(u"K", T))
        Λ_val = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(m, T))
        V_val = ustrip(u"Å^3", V)

        # Solve zV = (exp(βμ)/Λ³) · V for μ given a target zV.
        μ_for_zV(zV_target) = (log(zV_target * Λ_val^3 / V_val) / β) * u"eV"

        zV_targets = (0.5, 1.5, 3.0)
        μ_grid = [μ_for_zV(z) for z in zV_targets]

        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test)

        for (k, zV) in enumerate(zV_targets)
            @test isapprox(out.Xi[k, 1], exp(zV), rtol=1e-3)
            @test isapprox(out.mean_N[k, 1], zV, rtol=1e-3)
            @test isapprox(out.var_N[k, 1], zV, rtol=1e-3)
            @test isapprox(out.mean_U[k, 1], 0.0, atol=1e-12)
        end
    end

    @testset "Ξ monotone in μ at fixed T" begin
        N_max = 20
        N_values = collect(0:N_max)
        K = 120
        n_iters = 5000
        ω0_test = (K + 1) / K

        ns_outputs = [DataFrame(iter=collect(1:n_iters), emax=zeros(n_iters))
                      for _ in N_values]

        V = 1000.0u"Å^3"
        m = 40.0u"u"
        T_grid = [300.0u"K", 600.0u"K"]

        μ_grid = collect(range(-0.30u"eV", -0.20u"eV", length=10))

        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test)

        @test all(diff(out.Xi, dims=1) .> 0)
    end

    @testset "argument validation" begin
        K = 120
        n_iters = 100
        ns_outputs = [DataFrame(iter=collect(1:n_iters), emax=zeros(n_iters))
                      for _ in 1:3]

        # Missing N=0 in N_values
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            ns_outputs, [1, 2, 3],
            1000.0u"Å^3", 40.0u"u",
            [-0.25u"eV"], [300.0u"K"];
            n_walkers=K)

        # Length mismatch between ns_outputs and N_values
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            ns_outputs, [0, 1, 2, 3],
            1000.0u"Å^3", 40.0u"u",
            [-0.25u"eV"], [300.0u"K"];
            n_walkers=K)

        # Length mismatch between live_emax and N_values
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            ns_outputs, [0, 1, 2],
            1000.0u"Å^3", 40.0u"u",
            [-0.25u"eV"], [300.0u"K"];
            n_walkers=K,
            live_emax=[Float64[], Float64[]])
    end

end


@testset "Atomistic GC-NS, fixed-N construction (3D LJ fluid)" begin
    using Random

    # ===== System ==========================================================
    # 12 Å cubic periodic box with an Ar-like LJ fluid: ε = 0.01 eV (T_crit ~ 1.3ε/k_B
    # ≈ 150 K), σ = 2.5 Å, cutoff 3σ with energy shift. At T = 200 K we have
    # k_B T ≈ 1.72 ε — supercritical, safely away from coexistence physics.
    L = 12.0
    box = [[L * u"Å", 0u"Å", 0u"Å"],
           [0u"Å", L * u"Å", 0u"Å"],
           [0u"Å", 0u"Å", L * u"Å"]]
    V_box = L^3 * u"Å^3"
    mass = 40.0u"u"
    lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=3.0, shift=true)

    T_val = 200.0u"K"
    T_grid = [T_val]
    kb = 8.617333262e-5
    β = 1.0 / (kb * ustrip(u"K", T_val))

    N_values = collect(0:4)
    K = 48
    mc_steps = 500
    n_ns_steps = 800
    n_equi = 2000
    n_sample = 8000

    # Build N atoms on a sparse cubic-corner grid (max 5 spots used; for N ≤ 4
    # this gives initial separations ≈ L/2 = 6 Å, well outside the LJ wall).
    function _place_n_atoms(N::Int)
        spots = [(0.25, 0.25, 0.25),
                 (0.75, 0.25, 0.25),
                 (0.25, 0.75, 0.25),
                 (0.75, 0.75, 0.25),
                 (0.25, 0.25, 0.75)]
        positions = Pair{Symbol, Vector{Float64}}[]
        for i in 1:N
            x, y, z = spots[i]
            x += 0.02 * (rand() - 0.5)
            y += 0.02 * (rand() - 0.5)
            z += 0.02 * (rand() - 0.5)
            push!(positions, :Ar => [x, y, z])
        end
        return positions
    end

    function _build_liveset_N(N::Int)
        walkers = [AtomWalker(FastSystem(periodic_system(_place_n_atoms(N), box, fractional=true)))
                   for _ in 1:K]
        return LJAtomWalkers(walkers, lj)
    end

    function _run_ns_at_N(N::Int, seed::Int)
        Random.seed!(seed)
        liveset = _build_liveset_N(N)
        ns_params = NestedSamplingParameters(
            mc_steps=mc_steps,
            initial_step_size=0.3,
            step_size=0.3,
            step_size_lo=0.01,
            step_size_up=2.0,
            accept_range=(0.25, 0.75),
            allowed_fail_count=1000,
            energy_perturbation=1e-12,
        )
        save_strategy = SaveEveryN(
            df_filename="_test_lj_ns_$(N).csv",
            wk_filename="_test_lj_ns_$(N).traj.extxyz",
            ls_filename="_test_lj_ns_$(N).ls.extxyz",
            n_traj=100_000,
            n_snap=100_000,
            n_info=100_000,
        )
        df, final_liveset, _ = nested_sampling(
            liveset, ns_params, n_ns_steps, MCRandomWalkClone(), save_strategy)
        live_emax_N = [ustrip(u"eV", w.energy) for w in final_liveset.walkers]
        rm("_test_lj_ns_$(N).csv", force=true)
        rm("_test_lj_ns_$(N).traj.extxyz", force=true)
        rm("_test_lj_ns_$(N).ls.extxyz", force=true)
        return df, live_emax_N
    end

    function _run_nvt_at_N(N::Int, seed::Int)
        Random.seed!(seed)
        atoms_list = _place_n_atoms(N)
        sys = FastSystem(periodic_system(atoms_list, box, fractional=true))
        walker = AtomWalker(sys)
        mc_params = MetropolisMCParameters(
            [ustrip(u"K", T_val)];
            equilibrium_steps=n_equi,
            sampling_steps=n_sample,
            step_size=0.3,
            step_size_lo=0.05,
            step_size_up=1.0,
            accept_range=(0.3, 0.7),
            random_seed=seed,
        )
        energies, _, _, _ = monte_carlo_sampling(MCRandomWalkMaxE(), walker, lj, mc_params)
        return energies[1]
    end

    # ===== Run NS at each N ≥ 2 ==========================================
    # N = 0 is handled by the function's special case.
    # N = 1 cannot make NS progress: a single atom in a periodic box has no
    # pair interactions, so every walker has E = 0 exactly and every MC walk
    # is rejected by the energy-ceiling check. Z_NS^{(1)} = 1 by definition,
    # so we synthesize the DataFrame directly. Using a long, all-zero-energy
    # df with default ω₀ yields the same `K/(K+1)` multiplicative scaling as
    # the real NS runs at N ≥ 2, keeping the Ξ assembly internally consistent.
    ns_outputs = Vector{DataFrame}(undef, length(N_values))
    live_emax_all = Vector{Vector{Float64}}(undef, length(N_values))
    ns_outputs[1] = DataFrame(iter=Int[], emax=Float64[])
    live_emax_all[1] = Float64[]

    for (idx, N) in enumerate(N_values)
        N == 0 && continue
        if N == 1
            n_synth = 5000
            ns_outputs[idx] = DataFrame(iter=collect(1:n_synth), emax=zeros(n_synth))
            live_emax_all[idx] = zeros(Float64, K)
            continue
        end
        df, live_emax_N = _run_ns_at_N(N, 1000 + N)
        ns_outputs[idx] = df
        live_emax_all[idx] = live_emax_N
    end

    @testset "NS run schema and basic shape" begin
        for (idx, N) in enumerate(N_values)
            (N == 0 || N == 1) && continue
            df = ns_outputs[idx]
            @test names(df) == ["iter", "emax"]
            @test nrow(df) >= 0.5 * n_ns_steps
            @test length(live_emax_all[idx]) == K
        end
    end

    # ===== Per-N ⟨U⟩: NS vs NVT cross-check =============================
    # Construct per-N ⟨U⟩_NS with live-set tail closure, compare to NVT mean.
    # internal_energy uses the same ω·exp(−βE) weighting as the function.
    function _U_with_tail(df::DataFrame, live_emax_N::Vector{Float64})
        ωi = ωᵢ(df.iter, K)
        n_iters = length(df.iter)
        tail_w = (K / (K + 1))^n_iters / K
        ωi_full = vcat(ωi, fill(tail_w, K))
        Es_full = vcat(collect(Float64, df.emax), live_emax_N)
        return internal_energy(β, ωi_full, Es_full)
    end

    function _U_no_tail(df::DataFrame)
        ωi = ωᵢ(df.iter, K)
        Es = collect(Float64, df.emax)
        return internal_energy(β, ωi, Es)
    end

    # Cache NVT cross-check results — used by two test blocks.
    U_NVT_by_N = Dict{Int, Float64}()
    for N in N_values
        N == 0 && continue
        U_NVT_by_N[N] = _run_nvt_at_N(N, 2000 + N)
    end

    @testset "NS-vs-NVT ⟨U⟩_N agreement at T = 200 K" begin
        # The tolerance is deliberately loose. NS at small K with single-atom
        # random walks systematically underestimates |⟨U⟩| for N ≥ 3 in this
        # supercritical-but-cold regime (k_B T ≈ 1.7 ε): once the live set
        # collapses near the well bottom, the empirical df.emax sequence is too
        # coarse in the thermally relevant E band (a few k_B T above the
        # minimum) to faithfully reproduce the canonical distribution. The
        # observed bias is ~30–50 % at N = 3, 4, scaling with the number of
        # interacting pairs. NVT MC at n_sample = 8000 is converged to
        # micro-eV precision, so it is the trustworthy reference here.
        #
        # The test's purpose is to catch normalization / unit-conversion / sign
        # bugs in the per-N evidence assembly (which would produce errors of
        # 100 % or more), not to certify NS sampling quality.
        for (idx, N) in enumerate(N_values)
            N == 0 && continue
            U_NS = _U_with_tail(ns_outputs[idx], live_emax_all[idx])
            U_NVT = U_NVT_by_N[N]
            @test isapprox(U_NS, U_NVT; rtol=0.60, atol=3e-3)
            # Sign check: NS and NVT must both be non-positive in the
            # attractive regime (this catches a wrong-sign error that the
            # loose isapprox would otherwise miss for U ≈ 0).
            @test U_NS <= 1e-4
        end
    end

    @testset "Live-set tail correction moves U toward NVT" begin
        # For at least one N, the tail-corrected ⟨U⟩ should be no further from
        # NVT than the uncorrected ⟨U⟩. Tested at N = 2 where the well is
        # deepest in relative terms and live walkers concentrate below the
        # latest recorded emax. Includes a small slack to absorb the case
        # where the two estimates are essentially indistinguishable.
        N_check = 2
        idx_check = findfirst(==(N_check), N_values)
        df = ns_outputs[idx_check]
        live_emax_N = live_emax_all[idx_check]
        U_with = _U_with_tail(df, live_emax_N)
        U_without = _U_no_tail(df)
        U_NVT = U_NVT_by_N[N_check]
        @test abs(U_with - U_NVT) <= abs(U_without - U_NVT) + 1e-6
    end

    @testset "Per-N ⟨U⟩ decreases (more negative) with N" begin
        # Attractive LJ regime → adding a particle adds (on average) negative
        # pair contributions to ⟨U⟩_N. Allow 1 meV slack per step for noise.
        U_N_vec = Float64[]
        push!(U_N_vec, 0.0)  # N = 0 by definition
        for (idx, N) in enumerate(N_values)
            N == 0 && continue
            push!(U_N_vec, _U_with_tail(ns_outputs[idx], live_emax_all[idx]))
        end
        # N=1 → N=2: pair energy appears, U drops.
        # N=2 → N=3 → N=4: each adds more pair contributions.
        for i in 2:length(U_N_vec)
            @test U_N_vec[i] <= U_N_vec[i-1] + 1e-3
        end
    end

    # ===== End-to-end Ξ(μ, T) via gc_thermodynamic_stats_fixed_N =========
    # μ values placed in the dilute-gas window: at zV ≈ 1, μ ≈ -0.205 eV
    # (from the ideal-gas zV relation with V = 1728 Å³, m = 40u, T = 200 K).
    # μ_grid below brackets ⟨N⟩ values around 1 and 2.
    μ_grid = [-0.215u"eV", -0.200u"eV"]

    out = gc_thermodynamic_stats_fixed_N(
        ns_outputs, N_values, V_box, mass, μ_grid, T_grid;
        n_walkers=K, live_emax=live_emax_all)

    @testset "End-to-end Ξ assembly: returns expected shape" begin
        @test size(out.Xi) == (length(μ_grid), length(T_grid))
        @test size(out.mean_N) == (length(μ_grid), length(T_grid))
        @test size(out.var_N) == (length(μ_grid), length(T_grid))
        @test size(out.mean_U) == (length(μ_grid), length(T_grid))
        @test all(out.Xi .> 0)
        @test all(out.mean_N .>= 0)
        @test all(out.var_N .>= 0)
    end

    @testset "⟨N⟩ in expected range for chosen μ window" begin
        @test 0.3 < out.mean_N[1, 1] < 2.5
        @test 0.5 < out.mean_N[2, 1] < 3.5
    end

    @testset "Ξ and ⟨N⟩ monotone in μ at fixed T" begin
        @test out.Xi[2, 1] > out.Xi[1, 1]
        @test out.mean_N[2, 1] > out.mean_N[1, 1]
    end

    @testset "Grand ⟨U⟩ is non-positive (attractive regime, low density)" begin
        # In the dilute attractive regime, ⟨U⟩_GC should be ≤ 0 (with slack
        # for sampling noise at the lowest μ where ⟨N⟩ ≈ 0).
        @test out.mean_U[1, 1] <= 1e-3
        @test out.mean_U[2, 1] <= 1e-3
    end
end
