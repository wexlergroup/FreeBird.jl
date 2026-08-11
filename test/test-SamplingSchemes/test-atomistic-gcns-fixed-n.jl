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

        # Duplicate entries in N_values
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            ns_outputs, [0, 1, 1],
            1000.0u"Å^3", 40.0u"u",
            [-0.25u"eV"], [300.0u"K"];
            n_walkers=K)

        # Non-positive temperature
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            ns_outputs, [0, 1, 2],
            1000.0u"Å^3", 40.0u"u",
            [-0.25u"eV"], [0.0u"K"];
            n_walkers=K)
    end

    @testset "empty_energy: offset invariance, empty-sector restoration, guards" begin
        # Deterministic synthetic ladders — no NS runs. Sector N carries a
        # linearly descending emax ladder ending at its floor, plus K live
        # energies at the floor.
        K = 64
        n_iters = 400
        ω0_test = (K + 1) / K
        N_values = collect(0:4)

        V = 1000.0u"Å^3"
        m = 40.0u"u"
        T_grid = [200.0u"K", 350.0u"K"]
        kb = 8.617333262e-5
        βs = [1.0 / (kb * ustrip(u"K", T)) for T in T_grid]
        Λs = [ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(m, T)) for T in T_grid]
        V_val = ustrip(u"Å^3", V)
        # μ targets are set at the first (coldest) temperature; the second
        # temperature exercises the per-T empty-sector row at the same μ.
        μ_for_zV(zV_target) = (log(zV_target * Λs[1]^3 / V_val) / βs[1]) * u"eV"

        E_floor(N) = -0.02 * N
        ladder(N) = E_floor(N) .+ 0.05 .* (1.0 .- collect(1:n_iters) ./ n_iters)
        empty_df() = DataFrame(iter=Int[], emax=Float64[])
        ns_outputs = [N == 0 ? empty_df() :
                      DataFrame(iter=collect(1:n_iters), emax=ladder(N))
                      for N in N_values]
        live_all = [N == 0 ? Float64[] : fill(E_floor(N), K) for N in N_values]

        μ_grid = [μ_for_zV(0.3), μ_for_zV(1.0)]

        out0 = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, live_emax=live_all)

        @testset "offset invariance" begin
            # Shifting every ledger and live energy by a constant E0 while
            # declaring the same E0 as the empty-sector energy leaves ⟨N⟩ and
            # Var N unchanged, scales Ξ by exp(−βE0), and offsets ⟨U⟩ by
            # exactly E0.
            E0 = -1.5
            ns_shift = [N == 0 ? empty_df() :
                        DataFrame(iter=collect(1:n_iters), emax=ladder(N) .+ E0)
                        for N in N_values]
            live_shift = [N == 0 ? Float64[] : fill(E_floor(N) + E0, K)
                          for N in N_values]

            outE = gc_thermodynamic_stats_fixed_N(
                ns_shift, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_shift, empty_energy=E0)

            @test all(abs.(outE.mean_N .- out0.mean_N) .<= 1e-12)
            @test all(abs.(outE.var_N .- out0.var_N) .<= 1e-12)
            for j in eachindex(T_grid)
                @test all(isapprox.(outE.Xi[:, j], out0.Xi[:, j] .* exp(-βs[j] * E0),
                                    rtol=1e-10))
            end
            @test all(abs.(outE.mean_U .- (out0.mean_U .+ E0)) .<= 1e-12)
        end

        @testset "empty-sector restoration" begin
            # Constant-energy sectors E ≡ E0 + N·u make the truncated grand sum
            # exactly computable: Z_NS^{(N)} = M · exp(−β(E0 + N·u)) for N ≥ 1,
            # with the prior mass M = Σᵢ ωᵢ + tail replicated analytically. The
            # default call reproduces the biased sum in which the empty sector
            # enters with weight 1 instead of exp(−βE0); empty_energy = E0
            # reproduces the exact sum.
            E0 = -1.5
            u1 = -0.013
            E_N(N) = E0 + N * u1
            ns_const = [N == 0 ? empty_df() :
                        DataFrame(iter=collect(1:n_iters), emax=fill(E_N(N), n_iters))
                        for N in N_values]
            live_const = [N == 0 ? Float64[] : fill(E_N(N), K) for N in N_values]

            r = K / (K + 1)
            M_dead = sum(ω0_test * (1 / (K + 1)) * r^i for i in 1:n_iters)
            M_tail = ω0_test * r^n_iters
            M = M_dead + M_tail

            out_def = gc_thermodynamic_stats_fixed_N(
                ns_const, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_const)
            out_fix = gc_thermodynamic_stats_fixed_N(
                ns_const, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_const, empty_energy=E0)

            wavg(w, x) = sum(w .* x) / sum(w)
            Ns = collect(0:4)

            for (k, μ) in enumerate(μ_grid), j in eachindex(T_grid)
                βj = βs[j]
                zV = exp(βj * ustrip(u"eV", μ)) * V_val / Λs[j]^3
                terms = [zV^N / factorial(N) * M * exp(-βj * E_N(N)) for N in 1:4]
                w_ex = vcat(exp(-βj * E0), terms)  # exact: N = 0 at exp(−βE0)
                w_b = vcat(1.0, terms)             # biased: N = 0 at weight 1
                Es_ex = vcat(E0, [E_N(N) for N in 1:4])
                Es_b = vcat(0.0, [E_N(N) for N in 1:4])  # default: mean_E(0) = 0

                @test isapprox(out_fix.Xi[k, j], sum(w_ex), rtol=1e-10)
                @test isapprox(out_fix.mean_N[k, j], wavg(w_ex, Ns), rtol=1e-10)
                @test isapprox(out_fix.var_N[k, j],
                               wavg(w_ex, Ns .^ 2) - wavg(w_ex, Ns)^2, rtol=1e-10)
                @test isapprox(out_fix.mean_U[k, j], wavg(w_ex, Es_ex), rtol=1e-10)

                @test isapprox(out_def.Xi[k, j], sum(w_b), rtol=1e-10)
                @test isapprox(out_def.mean_N[k, j], wavg(w_b, Ns), rtol=1e-10)
                @test isapprox(out_def.var_N[k, j],
                               wavg(w_b, Ns .^ 2) - wavg(w_b, Ns)^2, rtol=1e-10)
                @test isapprox(out_def.mean_U[k, j], wavg(w_b, Es_b), rtol=1e-10)

                # The bias is material in the cold dilute column: the default
                # overstates ⟨N⟩ (analytic gaps 0.713 and 0.277 at 200 K).
                if j == 1
                    @test out_def.mean_N[k, 1] > out_fix.mean_N[k, 1] + 0.05
                end
            end
        end

        @testset "guards" begin
            @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
                ns_outputs, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_all, empty_energy=NaN)
            @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
                ns_outputs, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_all, empty_energy=Inf)
            @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
                ns_outputs, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_all, empty_energy=-Inf)
        end
    end

end


@testset "Atomistic GC-NS, fixed-N parity fields" begin
    # Deterministic synthetic ladders — no NS runs. Every closure is gated
    # against an independent reference: direct sample-level reweighting of the
    # raw ledger weights accumulated in BigFloat (no energy shift needed at
    # that precision), a different decomposition than the library's
    # law-of-total-expectation assembly.
    K = 64
    n_iters = 400
    ω0_test = (K + 1) / K
    N_values = collect(0:4)

    V = 1000.0u"Å^3"
    m = 40.0u"u"
    T_grid = [200.0u"K", 350.0u"K"]
    kb = 8.617333262e-5
    βs = [1.0 / (kb * ustrip(u"K", T)) for T in T_grid]
    Λs = [ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(m, T)) for T in T_grid]
    V_val = ustrip(u"Å^3", V)
    μ_for_zV(zV_target) = (log(zV_target * Λs[1]^3 / V_val) / βs[1]) * u"eV"
    μ_grid = [μ_for_zV(0.3), μ_for_zV(1.0)]
    empty_df() = DataFrame(iter=Int[], emax=Float64[])
    r = K / (K + 1)

    ladder(base, N) = base + N * (-0.02) .+ 0.05 .* (1.0 .- collect(1:n_iters) ./ n_iters)
    floor_of(base, N) = base + N * (-0.02)

    function brute_reference(ns_outs, live_all, N_vals, empty_E, obs_cols,
                             live_obs, μ, T)
        β = 1.0 / (kb * ustrip(u"K", T))
        Λ = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(m, T))
        log_zV = β * ustrip(u"eV", μ) + log(V_val) - 3 * log(Λ)
        S0 = BigFloat(0); SN = BigFloat(0); SN2 = BigFloat(0)
        SE = BigFloat(0); SE2 = BigFloat(0); SEN = BigFloat(0)
        SA = Dict(col => BigFloat(0) for col in obs_cols)
        function add!(E, lw, N, avals)
            w = exp(BigFloat(lw) - BigFloat(β) * BigFloat(E) +
                    N * BigFloat(log_zV) - log(BigFloat(factorial(N))))
            S0 += w; SN += w * N; SN2 += w * N^2
            SE += w * E; SE2 += w * E^2; SEN += w * E * N
            for col in obs_cols
                SA[col] += w * avals[col]
            end
            return nothing
        end
        for (i, N) in enumerate(N_vals)
            if N == 0
                a0 = Dict(col => sum(Float64, live_obs[i][col]) / length(live_obs[i][col])
                          for col in obs_cols)
                add!(empty_E, 0.0, 0, a0)
                continue
            end
            df = ns_outs[i]
            for row in 1:nrow(df)
                lw = log(ω0_test) + log(1 / (K + 1)) + df.iter[row] * log(r)
                arow = Dict(col => Float64(df[row, col]) for col in obs_cols)
                add!(df.emax[row], lw, N, arow)
            end
            n_it = maximum(df.iter)
            lw_tail = log(ω0_test) + n_it * log(r) - log(K)
            for (s, E) in enumerate(live_all[i])
                atail = Dict(col => Float64(live_obs === nothing ? 0.0 :
                                            live_obs[i][col][s]) for col in obs_cols)
                add!(E, lw_tail, N, atail)
            end
        end
        mN = SN / S0; mE = SE / S0
        return (logXi=Float64(log(S0)), mean_N=Float64(mN),
                var_N=Float64(SN2 / S0 - mN^2), mean_U=Float64(mE),
                var_U=Float64(SE2 / S0 - mE^2),
                cov_UN=Float64(SEN / S0 - mE * mN),
                obs=Dict(col => Float64(SA[col] / S0) for col in obs_cols))
    end

    @testset "brute-force closure, ordinary and offset ladders" begin
        # Variant (a): energies of ordinary magnitude, default empty_energy.
        # Variant (b): every energy carries a -1000 eV offset with
        # empty_energy to match — the E_shift cancellation guard is the only
        # thing keeping var_U/cov_UN at rtol 1e-10 there.
        for (base, E0) in ((0.0, 0.0), (-1000.0, -1000.0))
            ns_outputs = [N == 0 ? empty_df() :
                          DataFrame(iter=collect(1:n_iters), emax=ladder(base, N))
                          for N in N_values]
            live_all = [N == 0 ? Float64[] : fill(floor_of(base, N), K)
                        for N in N_values]
            out = gc_thermodynamic_stats_fixed_N(
                ns_outputs, N_values, V, m, μ_grid, T_grid;
                n_walkers=K, ω0=ω0_test, live_emax=live_all, empty_energy=E0)
            # var_U/cov_UN at the 10³ eV offset sit on the assembly's
            # unshifted-first-moment floor (u_avg accumulates ~|E|·eps of
            # absolute error before the E_shift subtraction; the lattice
            # method shares the structure): measured ≤ 1.5e-8 relative here,
            # against ≈5e-6 with the cancellation guard disabled — the guard
            # stays load-bearing even at the relaxed gate. The offset variant
            # therefore gates var_U/cov_UN at 3e-7 and mean_N/var_N at 1e-9;
            # everything else holds 1e-10.
            rtol_uu = base == 0.0 ? 1e-10 : 3e-7
            rtol_nn = base == 0.0 ? 1e-10 : 1e-9
            for (k, μ) in enumerate(μ_grid), (j, T) in enumerate(T_grid)
                ref = brute_reference(ns_outputs, live_all, N_values, E0,
                                      Symbol[], nothing, μ, T)
                @test isapprox(out.logXi[k, j], ref.logXi, rtol=1e-10)
                @test isapprox(out.mean_N[k, j], ref.mean_N, rtol=rtol_nn)
                @test isapprox(out.var_N[k, j], ref.var_N, rtol=rtol_nn)
                @test isapprox(out.mean_U[k, j], ref.mean_U, rtol=1e-10)
                @test isapprox(out.var_U[k, j], ref.var_U, rtol=rtol_uu)
                @test isapprox(out.cov_UN[k, j], ref.cov_UN, rtol=rtol_uu)
            end
        end
    end

    @testset "overflow: logXi finite where Xi is Inf" begin
        ns_outputs = [N == 0 ? empty_df() :
                      DataFrame(iter=collect(1:n_iters), emax=zeros(n_iters))
                      for N in N_values]
        live_all = [N == 0 ? Float64[] : zeros(K) for N in N_values]
        μ_deep = [μ_for_zV(exp(180.0))]
        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_deep, T_grid[1:1];
            n_walkers=K, ω0=ω0_test, live_emax=live_all)
        @test out.Xi[1, 1] == Inf
        ref = brute_reference(ns_outputs, live_all, N_values, 0.0,
                              Symbol[], nothing, μ_deep[1], T_grid[1])
        @test isapprox(out.logXi[1, 1], ref.logXi, rtol=1e-12)
    end

    @testset "P(N | μ, T) recipe from log_Z_N" begin
        # E ≡ 0 ladders: the truncated weights are (zV)^N/N! · M with the
        # prior mass M = Σᵢ ωᵢ + tail, and M = 1 for the empty sector.
        ns_outputs = [N == 0 ? empty_df() :
                      DataFrame(iter=collect(1:n_iters), emax=zeros(n_iters))
                      for N in N_values]
        live_all = [N == 0 ? Float64[] : zeros(K) for N in N_values]
        M_dead = sum(ω0_test * (1 / (K + 1)) * r^i for i in 1:n_iters)
        M_tail = ω0_test * r^n_iters
        M = M_dead + M_tail
        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, live_emax=live_all)
        @test out.N_values == N_values
        @test size(out.log_Z_N) == (length(N_values), length(T_grid))
        for (k, μ) in enumerate(μ_grid), (j, T) in enumerate(T_grid)
            β = βs[j]
            logP = out.log_Z_N[:, j] .+ β .* ustrip(u"eV", μ) .* out.N_values
            P = exp.(logP .- maximum(logP))
            P ./= sum(P)
            zV = exp(β * ustrip(u"eV", μ)) * V_val / Λs[j]^3
            w_ref = [N == 0 ? 1.0 : zV^N / factorial(N) * M for N in N_values]
            @test isapprox(P, w_ref ./ sum(w_ref), rtol=1e-12)
            @test isapprox(sum(P .* out.N_values), out.mean_N[k, j], rtol=1e-12)
            @test isapprox(sum(P .* out.N_values .^ 2) - sum(P .* out.N_values)^2,
                           out.var_N[k, j], rtol=1e-12)
        end
    end

    @testset "observables passthrough" begin
        obs_cols = [:A, :B]
        ns_outputs = Vector{DataFrame}(undef, length(N_values))
        live_all = Vector{Vector{Float64}}(undef, length(N_values))
        live_obs = Vector{Dict{Symbol,Vector{Float64}}}(undef, length(N_values))
        for (i, N) in enumerate(N_values)
            if N == 0
                ns_outputs[i] = empty_df()
                live_all[i] = Float64[]
                live_obs[i] = Dict(:A => [7.0, 8.0], :B => [0.0])
                continue
            end
            Es = ladder(0.0, N)
            ns_outputs[i] = DataFrame(iter=collect(1:n_iters), emax=Es,
                                      A=2.0 .* Es .+ 1.0, B=fill(Float64(N), n_iters))
            live_all[i] = fill(floor_of(0.0, N), K)
            live_obs[i] = Dict(:A => 2.0 .* live_all[i] .+ 1.0,
                               :B => fill(Float64(N), K))
        end
        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, live_emax=live_all,
            observable_cols=obs_cols, live_observables=live_obs)
        for (k, μ) in enumerate(μ_grid), (j, T) in enumerate(T_grid)
            ref = brute_reference(ns_outputs, live_all, N_values, 0.0,
                                  obs_cols, live_obs, μ, T)
            @test isapprox(out.observables[:A][k, j], ref.obs[:A], rtol=1e-10)
            @test isapprox(out.observables[:B][k, j], ref.obs[:B], rtol=1e-10)
            # B ≡ N per sector, except the empty sector's declared 0.0 — the
            # grand average must therefore reproduce ⟨N⟩ exactly.
            @test isapprox(out.observables[:B][k, j], out.mean_N[k, j], rtol=1e-12)
        end
        # Guards, mirroring the lattice method's validation.
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, live_emax=live_all, live_observables=live_obs)
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, observable_cols=obs_cols)
        bad_obs = [Dict(:A => d[:A]) for d in live_obs]
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, live_emax=live_all,
            observable_cols=obs_cols, live_observables=bad_obs)
    end

    @testset "cross-method pin on the coinciding two-sector ladder" begin
        # N_values = [0, 1] with E ≡ 0 is the one ladder where the atomistic
        # (zV/Λ³)^N/N! and lattice C(M,N)·e^{βμN} priors coincide exactly
        # under zV = M·e^{βμ_lat}.
        M_sites = 8
        N2 = [0, 1]
        ns2 = [empty_df(), DataFrame(iter=collect(1:n_iters), emax=zeros(n_iters))]
        live2 = [Float64[], zeros(K)]
        T1 = T_grid[1:1]
        μ_lat = [-0.05u"eV"]
        μ_atom = [μ_lat[1] + (log(M_sites * Λs[1]^3 / V_val) / βs[1]) * u"eV"]
        out_lat = gc_thermodynamic_stats_fixed_N(
            ns2, N2, M_sites, μ_lat, T1;
            n_walkers=K, ω0=ω0_test, live_emax=live2)
        out_atom = gc_thermodynamic_stats_fixed_N(
            ns2, N2, V, m, μ_atom, T1;
            n_walkers=K, ω0=ω0_test, live_emax=live2)
        @test abs(out_atom.mean_N[1, 1] - out_lat.mean_N[1, 1]) <= 1e-12
        @test abs(out_atom.var_N[1, 1] - out_lat.var_N[1, 1]) <= 1e-12
        @test isapprox(out_atom.logXi[1, 1], out_lat.logXi[1, 1], rtol=1e-12)
    end

    @testset "leading positions and legacy expressions" begin
        ns_outputs = [N == 0 ? empty_df() :
                      DataFrame(iter=collect(1:n_iters), emax=ladder(0.0, N))
                      for N in N_values]
        live_all = [N == 0 ? Float64[] : fill(floor_of(0.0, N), K)
                    for N in N_values]
        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V, m, μ_grid, T_grid;
            n_walkers=K, ω0=ω0_test, live_emax=live_all)
        @test propertynames(out) == (:Xi, :mean_N, :var_N, :mean_U, :logXi,
                                     :var_U, :cov_UN, :log_Z_N, :N_values,
                                     :observables)
        @test out[1] === out.Xi && out[2] === out.mean_N &&
              out[3] === out.var_N && out[4] === out.mean_U
        @test isempty(out.observables)
        # The four legacy fields must be bit-identical to the pre-extension
        # assembly, replicated here expression for expression.
        n = length(N_values)
        log_Z_NS = Matrix{Float64}(undef, n, length(T_grid))
        mean_E_N = Matrix{Float64}(undef, n, length(T_grid))
        for (i, N) in enumerate(N_values)
            if N == 0
                log_Z_NS[i, :] .= 0.0
                mean_E_N[i, :] .= 0.0
                continue
            end
            log_Z_NS[i, :], mean_E_N[i, :] = FreeBird.AnalysisTools._fixed_N_log_evidence(
                ns_outputs[i], T_grid;
                n_walkers=K, n_cull=1, ω0=ω0_test,
                live_energies=live_all[i], kb=kb)
        end
        log_fact = [FreeBird.AnalysisTools._log_factorial(N) for N in N_values]
        for (j, T) in enumerate(T_grid)
            β = 1.0 / (kb * ustrip(u"K", T))
            Λ_val = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(m, T))
            log_V_over_Λ3 = log(V_val) - 3 * log(Λ_val)
            for (k, μ) in enumerate(μ_grid)
                log_zV = β * ustrip(u"eV", μ) + log_V_over_Λ3
                log_w = log_Z_NS[:, j] .+ N_values .* log_zV .- log_fact
                max_log = maximum(log_w)
                ws = exp.(log_w .- max_log)
                sum_w = sum(ws)
                @test out.Xi[k, j] === exp(max_log) * sum_w
                @test out.mean_N[k, j] === sum(ws .* N_values) / sum_w
                mean_N2 = sum(ws .* (N_values .^ 2)) / sum_w
                @test out.var_N[k, j] === mean_N2 - (sum(ws .* N_values) / sum_w)^2
                @test out.mean_U[k, j] === sum(ws .* view(mean_E_N, :, j)) / sum_w
            end
        end
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

    # Uniform-prior NS initialization (i.i.d. random placement, overlaps
    # rejected), as canonical NS requires. A fixed-site initial live set is not
    # a prior draw and biases the NS evidence; see scripts/fixed_n_init_test.jl.
    function _uniform_walker_N(N::Int)
        while true
            coor = [:Ar => [rand(), rand(), rand()] for _ in 1:N]
            sys = FastSystem(periodic_system(coor, box, fractional=true))
            E = ustrip(u"eV", interacting_energy(sys, lj))
            (isfinite(E) && E < 100.0) && return AtomWalker(sys)
        end
    end

    function _build_liveset_N(N::Int)
        walkers = [_uniform_walker_N(N) for _ in 1:K]
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
            @test names(df) == ["iter", "emax", "log_compression"]
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
        # NS live sets are drawn from the uniform prior (i.i.d. random
        # placement), as canonical NS requires; the earlier fixed-site init
        # biased ⟨U⟩ low by ~40–50 %. With correct init the residual NS-vs-NVT
        # difference is sampling noise at these deliberately cheap CI params
        # (K = 48, n_ns = 800, one seed per N): deeper runs (K = 192,
        # n_ns = 8000) agree with NVT and an independent importance-sampling
        # anchor to ~0.1 % (see scripts/fixed_n_init_test.jl). NVT MC at
        # n_sample = 8000 is the trustworthy reference here.
        #
        # The bound below absorbs that cheap-params sampling noise (observed
        # ≤ 26 % here) while still catching normalization / unit / sign bugs in
        # the per-N evidence assembly (which produce > 100 % errors).
        for (idx, N) in enumerate(N_values)
            N == 0 && continue
            U_NS = _U_with_tail(ns_outputs[idx], live_emax_all[idx])
            U_NVT = U_NVT_by_N[N]
            @test isapprox(U_NS, U_NVT; rtol=0.40, atol=3e-3)
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


@testset "Atomistic GC-NS, fixed-N construction (LJ on LJ(111) surface)" begin
    using Random

    # ===== Substrate: LJ(111) slab =======================================
    # 2x2 in-plane rectangular supercell, 3 layers ABC stacking, all frozen.
    # In-plane NN distance a = 2^(1/6) σ (the LJ minimum). The rectangular
    # primitive cell on a (111) plane has dimensions (a, a√3) with 2 atoms
    # per cell; ABC stacking places adjacent layers at the up- and down-
    # triangle centroids of the layer below. Inter-layer spacing a√(2/3)
    # is the FCC (111) value. Box-x is technically shorter than 2·cutoff
    # (same convention as the 3D LJ fluid block above) so PBC self-image
    # contributions are present but appear as a constant offset that both
    # NS and NVT see identically — the cross-check remains valid.
    σ_val = 2.5
    ε_val = 0.01
    a = 2^(1/6) * σ_val
    Δz_layer = a * sqrt(2/3)
    Lx = 2.0 * a
    Ly = 2.0 * a * sqrt(3.0)
    Lz = 22.0  # vacuum above 3-layer slab > 2·cutoff = 15 Å

    # Layer-internal positions, given as fractions of the primitive
    # rectangle (Lx/2, Ly/2): A = layer with 2 atoms at primitive corners;
    # B and C are the up- and down-triangle centroids of A.
    layer_offsets = [
        [(0.0, 0.0), (0.5, 0.5)],         # A
        [(0.5, 1.0/6.0), (0.0, 2.0/3.0)], # B
        [(0.5, 5.0/6.0), (0.0, 1.0/3.0)], # C
    ]

    function _build_substrate_positions()
        positions = Pair{Symbol, Vector{Float64}}[]
        for (l, layer) in enumerate(layer_offsets)
            z_frac = (l - 1) * Δz_layer / Lz
            for (fx_prim, fy_prim) in layer
                for ix in 0:1, iy in 0:1
                    fx = (ix + fx_prim) / 2.0
                    fy = (iy + fy_prim) / 2.0
                    push!(positions, :Ar => [fx, fy, z_frac])
                end
            end
        end
        return positions
    end

    box = [[Lx * u"Å", 0u"Å", 0u"Å"],
           [0u"Å", Ly * u"Å", 0u"Å"],
           [0u"Å", 0u"Å", Lz * u"Å"]]
    V_box = (Lx * Ly * Lz) * u"Å^3"
    mass = 40.0u"u"

    lj = LJParameters(epsilon=ε_val, sigma=σ_val, cutoff=3.0, shift=true)
    # 2-component CPS for adsorbate × surface: [aa, as, ss] all identical.
    ljs = CompositeParameterSets(2, [lj, lj, lj])

    substrate_positions = _build_substrate_positions()
    substrate_sys = FastSystem(periodic_system(substrate_positions, box, fractional=true))
    surface = AtomWalker(substrate_sys; freeze_species=[:Ar])
    surface.energy_frozen_part = interacting_energy(surface.configuration, lj)
    E_surf_self = ustrip(u"eV", surface.energy_frozen_part)

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

    # ----- Adsorbate placement: above hollow-ish xy sites, well above the
    # top substrate layer (z = 2 Δz_layer + 1.5 σ).
    top_z_frac = 2 * Δz_layer / Lz
    init_z_frac = top_z_frac + (1.5 * σ_val) / Lz
    hollow_sites = [(0.25, 0.25),
                    (0.75, 0.25),
                    (0.25, 0.75),
                    (0.75, 0.75)]

    function _place_n_adsorbates(N::Int)
        positions = Pair{Symbol, Vector{Float64}}[]
        for i in 1:N
            fx, fy = hollow_sites[i]
            fx += 0.02 * (rand() - 0.5)
            fy += 0.02 * (rand() - 0.5)
            fz = init_z_frac + 0.01 * (rand() - 0.5)
            push!(positions, :Ar => [fx, fy, fz])
        end
        return positions
    end

    # Uniform-prior NS initialization over the full box (see the 3D LJ block).
    # Adsorbates landing in/near the slab are rejected here / discarded by NS.
    function _uniform_walker_N(N::Int)
        while true
            coor = [:Ar => [rand(), rand(), rand()] for _ in 1:N]
            w = AtomWalker(FastSystem(periodic_system(coor, box, fractional=true)))
            E_ads = ustrip(u"eV", interacting_energy(w.configuration, ljs,
                        w.list_num_par, w.frozen, surface.configuration))
            (isfinite(E_ads) && E_ads < 100.0) && return w
        end
    end

    function _build_liveset_N(N::Int)
        walkers = [_uniform_walker_N(N) for _ in 1:K]
        return LJSurfaceWalkers(walkers, ljs, surface)
    end

    function _run_ns_at_N_surface(N::Int, seed::Int)
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
            df_filename="_test_surf_ns_$(N).csv",
            wk_filename="_test_surf_ns_$(N).traj.extxyz",
            ls_filename="_test_surf_ns_$(N).ls.extxyz",
            n_traj=100_000,
            n_snap=100_000,
            n_info=100_000,
        )
        df, final_liveset, _ = nested_sampling(
            liveset, ns_params, n_ns_steps, MCRandomWalkClone(), save_strategy)
        live_emax_N = [ustrip(u"eV", w.energy) for w in final_liveset.walkers]
        rm("_test_surf_ns_$(N).csv", force=true)
        rm("_test_surf_ns_$(N).traj.extxyz", force=true)
        rm("_test_surf_ns_$(N).ls.extxyz", force=true)
        return df, live_emax_N
    end

    # Bespoke NVT MC loop with surface support. The production
    # `nvt_monte_carlo` does not yet take a surface argument; this mirrors
    # the surface-aware `MC_random_walk!` from `MonteCarloMoves` but
    # accepts on the Metropolis-Hastings ratio at temperature T.
    function _run_nvt_at_N_surface(N::Int, seed::Int)
        Random.seed!(seed)
        sys = FastSystem(periodic_system(_place_n_adsorbates(N), box, fractional=true))
        walker = AtomWalker(sys)
        walker.energy = interacting_energy(
            walker.configuration, ljs, walker.list_num_par, walker.frozen,
            surface.configuration) + surface.energy_frozen_part

        step_size = 0.3
        kbT_eV = kb * ustrip(u"K", T_val)

        function inner!(n_steps)
            energies = Vector{Float64}(undef, n_steps)
            for i in 1:n_steps
                config = walker.configuration
                free_index = FreeBird.MonteCarloMoves.free_par_index(walker)
                if isempty(free_index)
                    energies[i] = ustrip(u"eV", walker.energy)
                    continue
                end
                i_at = rand(free_index)
                prewalk = FreeBird.EnergyEval.single_site_energy(
                    i_at, config, ljs, walker.list_num_par, surface.configuration)
                pos = position(config, i_at)
                orig = pos
                pos = FreeBird.MonteCarloMoves.single_atom_random_walk!(pos, step_size)
                pos = FreeBird.MonteCarloMoves.periodic_boundary_wrap!(pos, config)
                config.position[i_at] = pos
                postwalk = FreeBird.EnergyEval.single_site_energy(
                    i_at, config, ljs, walker.list_num_par, surface.configuration)
                ΔE = postwalk - prewalk
                if ΔE < 0.0u"eV" || rand() < exp(-ustrip(u"eV", ΔE) / kbT_eV)
                    walker.energy = walker.energy + ΔE
                else
                    config.position[i_at] = orig
                end
                energies[i] = ustrip(u"eV", walker.energy)
            end
            return energies
        end

        inner!(n_equi)
        energies_sample = inner!(n_sample)
        return mean(energies_sample)
    end

    # ===== Run NS at each N >= 1 ========================================
    # N = 0 is handled by gc_thermodynamic_stats_fixed_N's special case.
    # N = 1 is non-trivial here (single adsorbate sees the substrate's
    # attractive potential), unlike the 3D LJ fluid block where N = 1 has
    # E ≡ 0 and NS cannot make progress.
    ns_outputs = Vector{DataFrame}(undef, length(N_values))
    live_emax_all = Vector{Vector{Float64}}(undef, length(N_values))
    ns_outputs[1] = DataFrame(iter=Int[], emax=Float64[])
    live_emax_all[1] = Float64[]

    for (idx, N) in enumerate(N_values)
        N == 0 && continue
        df, live_emax_N = _run_ns_at_N_surface(N, 3000 + N)
        ns_outputs[idx] = df
        live_emax_all[idx] = live_emax_N
    end

    @testset "NS run schema and basic shape" begin
        for (idx, N) in enumerate(N_values)
            N == 0 && continue
            df = ns_outputs[idx]
            @test names(df) == ["iter", "emax", "log_compression"]
            @test nrow(df) >= 0.5 * n_ns_steps
            @test length(live_emax_all[idx]) == K
        end
    end

    # ===== Per-N ⟨U⟩: NS vs NVT cross-check =============================
    # Compare on the adsorbate contribution ⟨U⟩ − E_surf_self, since the
    # constant substrate self-energy dominates the absolute scale and would
    # otherwise make rtol-based agreement trivially pass.
    function _U_with_tail(df::DataFrame, live_emax_N::Vector{Float64})
        ωi = ωᵢ(df.iter, K)
        n_iters = length(df.iter)
        tail_w = (K / (K + 1))^n_iters / K
        ωi_full = vcat(ωi, fill(tail_w, K))
        Es_full = vcat(collect(Float64, df.emax), live_emax_N)
        return internal_energy(β, ωi_full, Es_full)
    end

    U_NVT_by_N = Dict{Int, Float64}()
    for N in N_values
        N == 0 && continue
        U_NVT_by_N[N] = _run_nvt_at_N_surface(N, 4000 + N)
    end

    @testset "NS-vs-NVT ⟨U⟩_N agreement at T = 200 K" begin
        # NS live sets are drawn from the uniform prior (see the 3D LJ block);
        # with correct init the residual NS-vs-NVT difference is cheap-params
        # sampling noise (one seed per N at K = 48, n_ns = 800), not a
        # systematic bias. Comparison is on the adsorbate part of ⟨U⟩ (total
        # minus the constant substrate self-energy), since the constant offset
        # would otherwise dwarf any disagreement. The bound still catches
        # normalization / unit / sign bugs (which produce > 100 % errors).
        for (idx, N) in enumerate(N_values)
            N == 0 && continue
            U_NS_ads = _U_with_tail(ns_outputs[idx], live_emax_all[idx]) - E_surf_self
            U_NVT_ads = U_NVT_by_N[N] - E_surf_self
            @test isapprox(U_NS_ads, U_NVT_ads; rtol=0.40, atol=3e-3)
            # Adsorbate part should be non-positive in the attractive
            # regime (single Ar binds to the substrate at ~−3ε).
            @test U_NS_ads <= 1e-4
        end
    end

    @testset "Per-N ⟨U⟩ decreases (more negative) with N" begin
        # Each adsorbate adds attractive surface and (for N ≥ 2)
        # attractive lateral contributions. Allow 1 meV slack per step.
        U_N_vec = Float64[]
        push!(U_N_vec, 0.0)  # adsorbate contribution at N = 0 is 0
        for (idx, N) in enumerate(N_values)
            N == 0 && continue
            push!(U_N_vec, _U_with_tail(ns_outputs[idx], live_emax_all[idx]) - E_surf_self)
        end
        for i in 2:length(U_N_vec)
            @test U_N_vec[i] <= U_N_vec[i-1] + 1e-3
        end
    end

    # ===== End-to-end Ξ(μ, T) via gc_thermodynamic_stats_fixed_N ========
    # V passed is the simulation-box volume = NS prior volume per particle
    # (NOT a slab or adsorption-region volume). The substrate's
    # contribution to the configurational integral is fully captured by
    # the Boltzmann factor exp(−βU) in the per-N NS evidence; no
    # geometric V reduction is needed.
    #
    # μ window placed in the dilute-coverage regime. With ε = 0.01 eV
    # and a single-adsorbate binding of roughly −3 ε at hollow sites, we
    # expect ⟨N⟩ ~ 1–3 across the chosen μ pair.
    μ_grid = [-0.235u"eV", -0.215u"eV"]

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

    @testset "⟨N⟩ in sub-monolayer window for chosen μ" begin
        # Monolayer capacity ≈ 8 hollow sites on the 2x2 (111) supercell;
        # sub-monolayer means ⟨N⟩ well below 8.
        @test 0.1 < out.mean_N[1, 1] < 4.5
        @test 0.1 < out.mean_N[2, 1] < 4.5
    end

    @testset "Ξ and ⟨N⟩ monotone in μ at fixed T" begin
        @test out.Xi[2, 1] > out.Xi[1, 1]
        @test out.mean_N[2, 1] > out.mean_N[1, 1]
    end
end


@testset "Atomistic GC-NS, fixed-N ideal-adsorbate closure (LJ fcc(100) slab)" begin
    using Random

    # ===== Substrate: fcc(100) slab, all frozen ==========================
    # 3×3 in-plane square supercell with side 3·2^(1/6)·σ (the LJ-minimum
    # spacing), 3 layers in ABA stacking with the fcc(100) interlayer
    # spacing a/√2, box height 25 Å, fully periodic. Adsorbate-adsorbate
    # interactions are switched OFF (ε_aa = 0 in the composite parameter
    # set), so each adsorbate sees only the frozen slab field U_s(r): an
    # ideal gas in an external field, exactly solvable by quadrature. This
    # makes the block the one end-to-end surface fixture that closes
    # against exact references rather than shape and monotonicity checks
    # alone.
    σ_val = 2.5
    ε_val = 0.01
    a = 2^(1/6) * σ_val
    Lx = 3.0 * a
    Ly = 3.0 * a
    Lz = 25.0
    dz_layer = a / sqrt(2.0)
    z_base = 4.0

    slab_positions = Pair{Symbol,Vector{Float64}}[]
    for l in 0:2
        off = (l % 2 == 0) ? 0.0 : 0.5
        z = z_base + l * dz_layer
        for i in 0:2, j in 0:2
            push!(slab_positions, :Ar => [(i + off) / 3.0, (j + off) / 3.0, z / Lz])
        end
    end

    box = [[Lx * u"Å", 0u"Å", 0u"Å"],
           [0u"Å", Ly * u"Å", 0u"Å"],
           [0u"Å", 0u"Å", Lz * u"Å"]]
    V_box = (Lx * Ly * Lz) * u"Å^3"
    mass = 40.0u"u"

    lj = LJParameters(epsilon=ε_val, sigma=σ_val, cutoff=3.0, shift=true)
    # [aa, as, ss] upper-triangle order: the aa channel gets ε = 0.
    lj_zero = LJParameters(epsilon=0.0, sigma=σ_val, cutoff=3.0, shift=true)
    ljs = CompositeParameterSets(2, [lj_zero, lj, lj])

    substrate_sys = FastSystem(periodic_system(slab_positions, box, fractional=true))
    surface = AtomWalker(substrate_sys; freeze_species=[:Ar])
    surface.energy_frozen_part = interacting_energy(surface.configuration, lj)
    E_slab = ustrip(u"eV", surface.energy_frozen_part)

    T_val = 200.0u"K"
    T_grid = [T_val]
    kb = 8.617333262e-5
    β = 1.0 / (kb * ustrip(u"K", T_val))

    # ===== Quadrature reference ==========================================
    # Midpoint rule over the box for the per-particle field integrals
    # I = ∫ e^{−βU_s} d³r, I_u = ∫ U_s e^{−βU_s} d³r, I_u2 = ∫ U_s² e^{−βU_s} d³r,
    # plus the prior-retention fraction m₁ = |{r : U_s(r) < E_thresh}| / V
    # (see the walker-init note below). The pair kernel replicates the
    # shifted, 3σ-truncated lj_energy exactly; U is clamped at 50 eV, where
    # e^{−βU} is already 0, so the U²e^{−βU} accumuland stays finite inside
    # the hard cores. Grid (32, 32, 96) is the smallest rung of the halving
    # ladder (16, 16, 48) → (32, 32, 96) whose refinement step changes I by
    # < 1e-3 relative (measured: 2.3×10⁻⁵ on I and < 1e-3 on both U-moments
    # against (40, 40, 120)); the coarser rung (16, 16, 48) still moves the
    # U-moments by several ×10⁻³.
    slab_xyz = Matrix{Float64}(undef, 3, length(substrate_sys))
    for i in 1:length(substrate_sys)
        slab_xyz[:, i] .= ustrip.(u"Å", position(substrate_sys, i))
    end

    function _slab_field_moments(slab_xyz::Matrix{Float64}, Lx, Ly, Lz, σ, ε, β,
                                 nx::Int, ny::Int, nz::Int, E_thresh::Float64)
        rc = 3.0 * σ
        shift_e = 4.0 * ε * ((σ / rc)^12 - (σ / rc)^6)
        mind(d, L) = min(abs(d), L - abs(d))
        npos = size(slab_xyz, 2)
        hx, hy, hz = Lx / nx, Ly / ny, Lz / nz
        dV = hx * hy * hz
        I = 0.0; Iu = 0.0; Iu2 = 0.0; n_below = 0
        for iz in 1:nz
            z = (iz - 0.5) * hz
            for ix in 1:nx
                x = (ix - 0.5) * hx
                for iy in 1:ny
                    y = (iy - 0.5) * hy
                    U = 0.0
                    @inbounds for s in 1:npos
                        dx = mind(x - slab_xyz[1, s], Lx)
                        dy = mind(y - slab_xyz[2, s], Ly)
                        dz = mind(z - slab_xyz[3, s], Lz)
                        r = sqrt(dx * dx + dy * dy + dz * dz)
                        if r <= rc
                            sr6 = (σ / r)^6
                            U += 4.0 * ε * (sr6 * sr6 - sr6) - shift_e
                        end
                    end
                    U < E_thresh && (n_below += 1)
                    Uc = U > 50.0 ? 50.0 : U
                    e = exp(-β * Uc)
                    I += e
                    Iu += Uc * e
                    Iu2 += Uc * Uc * e
                end
            end
        end
        return I * dV, Iu * dV, Iu2 * dV, n_below / (nx * ny * nz)
    end

    E_init_thresh = 1.0e9
    I_ref, Iu, Iu2, m1 = _slab_field_moments(slab_xyz, Lx, Ly, Lz, σ_val, ε_val, β,
                                             32, 32, 96, E_init_thresh)
    ubar = Iu / I_ref          # per-particle ⟨U_s⟩ under e^{−βU_s}/I
    sig_u2 = Iu2 / I_ref - ubar^2
    Λ = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass, T_val))
    kT = kb * ustrip(u"K", T_val)
    # μ placed so the occupancy parameter q = e^{βμ}·I/Λ³ — the Poisson
    # rate of the untruncated ideal-adsorbate grand ensemble — lands at
    # 0.5 (dilute) and 1.0; both fall near −0.22 eV here.
    μ_grid = [kT * log(0.5 * Λ^3 / I_ref), kT * log(1.0 * Λ^3 / I_ref)] .* u"eV"

    # ===== Canonical NS ladders, N = 1..4 ================================
    # Walker-init rejection threshold: 1e9 eV, deliberately NOT the 100 eV
    # used by the LJ(111) block above. Rejecting E_ads ≥ 100 eV draws
    # truncates the uniform prior by a per-particle fraction worth ≈ 0.153
    # nats, so every per-N evidence would sit ≈ 0.153·N nats above the
    # reference — the bias alone exceeds the N = 4 gate, and combined with
    # the NS noise it erodes the others. At 1e9 eV only the hard-core
    # region is excluded, and the retained fraction m₁ ≈ 0.9972 (measured
    # on the quadrature grid) is folded into the reference as −N·ln m₁.
    K = 64
    mc_steps = 400
    n_ns_steps = 1200
    N_values = collect(0:4)
    seedbase = 3000   # gates calibrated at seedbases 3000 and 4000; 3000 ships

    function _uniform_walker_ideal(N::Int)
        while true
            coor = [:Ar => [rand(), rand(), rand()] for _ in 1:N]
            w = AtomWalker(FastSystem(periodic_system(coor, box, fractional=true)))
            E_ads = ustrip(u"eV", interacting_energy(w.configuration, ljs,
                        w.list_num_par, w.frozen, surface.configuration))
            (isfinite(E_ads) && E_ads < E_init_thresh) && return w
        end
    end

    tmpdir = mktempdir()
    function _run_ns_ideal(N::Int, seed::Int)
        Random.seed!(seed)
        walkers = [_uniform_walker_ideal(N) for _ in 1:K]
        liveset = LJSurfaceWalkers(walkers, ljs, surface)
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
            df_filename=joinpath(tmpdir, "_test_ideal_ns_$(N).csv"),
            wk_filename=joinpath(tmpdir, "_test_ideal_ns_$(N).traj.extxyz"),
            ls_filename=joinpath(tmpdir, "_test_ideal_ns_$(N).ls.extxyz"),
            n_traj=100_000,
            n_snap=100_000,
            n_info=100_000,
        )
        df, final_liveset, _ = nested_sampling(
            liveset, ns_params, n_ns_steps, MCRandomWalkClone(), save_strategy)
        live_emax_N = [ustrip(u"eV", w.energy) for w in final_liveset.walkers]
        return df, live_emax_N
    end

    ns_outputs = Vector{DataFrame}(undef, length(N_values))
    live_emax_all = Vector{Vector{Float64}}(undef, length(N_values))
    ns_outputs[1] = DataFrame(iter=Int[], emax=Float64[])
    live_emax_all[1] = Float64[]
    for (idx, N) in enumerate(N_values)
        N == 0 && continue
        df, live_emax_N = _run_ns_ideal(N, seedbase + N)
        ns_outputs[idx] = df
        live_emax_all[idx] = live_emax_N
    end
    rm(tmpdir; recursive=true, force=true)

    # ===== Stitch and exact references ===================================
    # empty_energy = E_slab: the N = 0 sector is the bare slab at its
    # frozen self-energy, which puts all five sectors on one energy scale.
    out_fix = gc_thermodynamic_stats_fixed_N(
        ns_outputs, N_values, V_box, mass, μ_grid, T_grid;
        n_walkers=K, live_emax=live_emax_all, empty_energy=E_slab)

    # Truncated-Poisson reference over the SAME N = 0..4 sector set, with
    # the corrected rate q_corr = q/m₁ (the NS prior per particle is the
    # m₁ fraction of the box). Ideal adsorbates make the sector energy
    # E_slab + Σᵢuᵢ with uᵢ i.i.d. (mean ubar, variance sig_u2), so
    #   ⟨U⟩ = E_slab + ubar·⟨N⟩,
    #   var_U = ubar²·Var N + sig_u2·⟨N⟩   (law of total variance),
    #   cov_UN = ubar·Var N.
    q_phys = [exp(β * ustrip(u"eV", μ)) * I_ref / Λ^3 for μ in μ_grid]
    q_corr = q_phys ./ m1
    function _trunc_poisson(q::Float64)
        Ns = 0:4
        w = [q^N / factorial(N) for N in Ns]
        S = sum(w)
        p = w ./ S
        mN = sum(p .* Ns)
        vN = sum(p .* Ns .^ 2) - mN^2
        return (S=S, mean_N=mN, var_N=vN,
                mean_U=E_slab + ubar * mN,
                var_U=ubar^2 * vN + sig_u2 * mN,
                cov_UN=ubar * vN)
    end

    @testset "per-N evidence vs quadrature reference" begin
        # Exact: log Z_N = −βE_slab + N·ln(I/(Λ³m₁)) − ln N! (the log_Z_N
        # return already carries the N·ln(V/Λ³) − ln N! part, so the
        # reference is stated on the same scale). Gates: ≥ 3× the larger
        # two-seedbase |deviation| per N — measured 0.0112/0.1030 (N = 1),
        # 0.1928/0.1310 (N = 2), 0.2889/0.0120 (N = 3), 0.1551/0.0077
        # (N = 4) nats at seedbases 3000/4000, consistent with the √(3N/K)
        # NS noise scale.
        tol_logZ = [0.35, 0.60, 0.90, 0.50]
        for (i, N) in enumerate(N_values)
            N == 0 && continue
            exact = -β * E_slab + N * (log(I_ref / Λ^3) - log(m1)) -
                    sum(log, 1:N; init=0.0)
            @test abs(out_fix.log_Z_N[i, 1] - exact) <= tol_logZ[N]
        end
    end

    @testset "stitched grand stats vs truncated Poisson" begin
        # Gates: ≥ 3× the max two-seedbase |deviation| per stat across both
        # μ — measured maxima: logXi 0.0647 nats, mean_N 0.0850 rel,
        # var_N 0.1055 rel, mean_U 1.55e-3 eV, var_U 0.1490 rel,
        # cov_UN 0.1390 rel. mean_U is gated absolutely: E_slab dominates
        # its magnitude, so a relative gate on the total would be
        # insensitive to the adsorbate part.
        for k in eachindex(μ_grid)
            tp = _trunc_poisson(q_corr[k])
            @test abs(out_fix.logXi[k, 1] - (-β * E_slab + log(tp.S))) <= 0.20
            @test abs(out_fix.mean_N[k, 1] / tp.mean_N - 1) <= 0.26
            @test abs(out_fix.var_N[k, 1] / tp.var_N - 1) <= 0.32
            @test abs(out_fix.mean_U[k, 1] - tp.mean_U) <= 5.0e-3
            @test abs(out_fix.var_U[k, 1] / tp.var_U - 1) <= 0.45
            @test abs(out_fix.cov_UN[k, 1] / tp.cov_UN - 1) <= 0.42
        end
    end

    @testset "default empty sector overstates ⟨N⟩" begin
        # Without empty_energy the N = 0 sector enters the grand sum at
        # weight 1 instead of e^{−βE_slab} = e^{+85.9}: it effectively
        # vanishes, and ⟨N⟩ collapses toward the N ≥ 1 conditional mean.
        # Measured gaps above the corrected ⟨N⟩ at the dilute μ: +0.778 and
        # +0.810 at seedbases 3000/4000 (the corrected values are ≈ 0.54
        # and ≈ 0.50); the margin ships at 0.25, ≥ 3× inside the measured
        # gaps while still asserting a ≥ +45% bias on the corrected ⟨N⟩.
        out_def = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V_box, mass, μ_grid, T_grid;
            n_walkers=K, live_emax=live_emax_all)
        @test out_def.mean_N[1, 1] > out_fix.mean_N[1, 1] + 0.25
    end
end


@testset "Atomistic observables ledger and passthrough (LJ on LJ(111) surface)" begin
    using Random

    # Same LJ(111) slab geometry as the block above, at reduced cost
    # (N = 2 free adsorbates, K = 32, 300 iterations): every assertion here
    # is a deterministic recording or passthrough closure given the seeded
    # run, not a sampling-quality gate. Covered: the per-dead-point
    # `observables` hook on an atomistic canonical run, the dead-point
    # callback's bit-exact pairing with the ledger, and `observable_cols`
    # passthrough of the recorded column through the fixed-N assembly.
    σ_val = 2.5
    ε_val = 0.01
    a = 2^(1/6) * σ_val
    Δz_layer = a * sqrt(2/3)
    Lx = 2.0 * a
    Ly = 2.0 * a * sqrt(3.0)
    Lz = 22.0

    layer_offsets = [
        [(0.0, 0.0), (0.5, 0.5)],         # A
        [(0.5, 1.0/6.0), (0.0, 2.0/3.0)], # B
        [(0.5, 5.0/6.0), (0.0, 1.0/3.0)], # C
    ]

    substrate_positions = Pair{Symbol,Vector{Float64}}[]
    for (l, layer) in enumerate(layer_offsets)
        z_frac = (l - 1) * Δz_layer / Lz
        for (fx_prim, fy_prim) in layer
            for ix in 0:1, iy in 0:1
                push!(substrate_positions,
                      :Ar => [(ix + fx_prim) / 2.0, (iy + fy_prim) / 2.0, z_frac])
            end
        end
    end

    box = [[Lx * u"Å", 0u"Å", 0u"Å"],
           [0u"Å", Ly * u"Å", 0u"Å"],
           [0u"Å", 0u"Å", Lz * u"Å"]]
    V_box = (Lx * Ly * Lz) * u"Å^3"
    mass = 40.0u"u"

    lj = LJParameters(epsilon=ε_val, sigma=σ_val, cutoff=3.0, shift=true)
    ljs = CompositeParameterSets(2, [lj, lj, lj])

    substrate_sys = FastSystem(periodic_system(substrate_positions, box, fractional=true))
    surface = AtomWalker(substrate_sys; freeze_species=[:Ar])
    surface.energy_frozen_part = interacting_energy(surface.configuration, lj)

    T_val = 200.0u"K"
    kb = 8.617333262e-5
    β = 1.0 / (kb * ustrip(u"K", T_val))

    N_free = 2
    K = 32
    mc_steps = 500
    n_ns_steps = 300

    function _uniform_walker_obs(N::Int)
        while true
            coor = [:Ar => [rand(), rand(), rand()] for _ in 1:N]
            w = AtomWalker(FastSystem(periodic_system(coor, box, fractional=true)))
            E_ads = ustrip(u"eV", interacting_energy(w.configuration, ljs,
                        w.list_num_par, w.frozen, surface.configuration))
            (isfinite(E_ads) && E_ads < 100.0) && return w
        end
    end

    # z_com: configuration-mean z of the free adsorbates in Å as a Float64 —
    # the adsorption-height observable natural to a slab geometry.
    z_com = function (cfg)
        s = 0.0
        for i in 1:length(cfg)
            s += ustrip(u"Å", position(cfg, i)[3])
        end
        return s / length(cfg)
    end

    callback_log = Tuple{Int,Float64}[]
    record_dead = (iter, w) -> push!(callback_log, (iter, ustrip(u"eV", w.energy)))

    Random.seed!(42)
    walkers = [_uniform_walker_obs(N_free) for _ in 1:K]
    liveset = LJSurfaceWalkers(walkers, ljs, surface)
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
    tmpdir = mktempdir()
    save_strategy = SaveEveryN(
        df_filename=joinpath(tmpdir, "_test_obs_ns.csv"),
        wk_filename=joinpath(tmpdir, "_test_obs_ns.traj.extxyz"),
        ls_filename=joinpath(tmpdir, "_test_obs_ns.ls.extxyz"),
        n_traj=100_000,
        n_snap=100_000,
        n_info=100_000,
    )
    df, final_liveset, _ = nested_sampling(
        liveset, ns_params, n_ns_steps, MCRandomWalkClone(), save_strategy;
        observables=[:z_com => z_com],
        dead_point_callback=record_dead)
    rm(tmpdir; recursive=true, force=true)

    @testset "observables ledger column" begin
        @test "z_com" in names(df)
        @test eltype(df.z_com) === Float64
        @test all(isfinite, df.z_com)
        @test all(0.0 .< df.z_com .< Lz)
    end

    @testset "dead-point callback pairing" begin
        # One invocation per recorded dead point, in ledger order, with the
        # culled walker's energy matching the emax column bit-exactly.
        @test length(callback_log) == nrow(df)
        @test all(callback_log[i][1] == df.iter[i] for i in eachindex(callback_log))
        @test all(callback_log[i][2] === df.emax[i] for i in eachindex(callback_log))
    end

    @testset "observable passthrough vs direct reweighting" begin
        # Single-sector demonstration: the run's N = 2 sector plus the
        # empty sector (whose live_observables entry declares the empty
        # configuration's observable value, 0.0 here). The reference is a
        # direct log-space reweighting of the same ledger column — dead
        # points at ω_i·e^{−βE_i}, live tail at the shared tail weight —
        # accumulated in BigFloat: a different decomposition of the same
        # sum, so agreement is deterministic given the run.
        live_emax_N = [ustrip(u"eV", w.energy) for w in final_liveset.walkers]
        live_z = [z_com(w.configuration) for w in final_liveset.walkers]

        N_values = [0, 2]
        μ_grid = [-0.22u"eV"]
        T_grid = [T_val]
        ns_outputs = [DataFrame(iter=Int[], emax=Float64[]), df]
        live_emax_all = [Float64[], live_emax_N]
        live_obs = [Dict(:z_com => [0.0]), Dict(:z_com => live_z)]

        out = gc_thermodynamic_stats_fixed_N(
            ns_outputs, N_values, V_box, mass, μ_grid, T_grid;
            n_walkers=K, live_emax=live_emax_all,
            observable_cols=[:z_com], live_observables=live_obs)

        Λ = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass, T_val))
        log_zV = β * ustrip(u"eV", μ_grid[1]) + log(ustrip(u"Å^3", V_box)) - 3 * log(Λ)
        log_r = log(BigFloat(K) / (K + 1))
        S0 = BigFloat(1)          # N = 0 sector: weight e^{−β·0}·(zV)⁰/0! = 1
        SA = BigFloat(0)          # its declared z_com is 0.0
        sector = 2 * BigFloat(log_zV) - log(BigFloat(2))   # N·ln zV − ln N!, N = 2
        for row in 1:nrow(df)
            lw = df.iter[row] * log_r - log(BigFloat(K + 1))
            w = exp(lw - BigFloat(β) * BigFloat(df.emax[row]) + sector)
            S0 += w
            SA += w * df.z_com[row]
        end
        n_it = maximum(df.iter)
        lw_tail = n_it * log_r - log(BigFloat(K))
        for (s, E) in enumerate(live_emax_N)
            w = exp(lw_tail - BigFloat(β) * BigFloat(E) + sector)
            S0 += w
            SA += w * live_z[s]
        end
        z_direct = Float64(SA / S0)

        @test isapprox(out.observables[:z_com][1, 1], z_direct, rtol=1e-10)
    end
end
