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
    end

end
