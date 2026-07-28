@testset "Lattice GC-NS, fixed-N construction" begin

    kb = 8.617333262e-5 # eV/K

    # ================================================================
    @testset "_log_binomial" begin
        lb = FreeBird.AnalysisTools._log_binomial
        @test lb(16, 0) == 0.0
        @test lb(16, 16) == 0.0
        @test isapprox(lb(16, 8), log(binomial(16, 8)), rtol=1e-12)
        # Beyond Int64 overflow of binomial(M, N)
        @test isapprox(lb(200, 100),
                       Float64(log(binomial(big(200), big(100)))), rtol=1e-10)
        @test_throws ArgumentError lb(4, 5)
        @test_throws ArgumentError lb(4, -1)
    end

    # ================================================================
    @testset "argument validation" begin
        μs = [0.0u"eV"]
        Ts = [300.0u"K"]
        empty_df() = DataFrame(iter=Int[], emax=Float64[])
        good_df() = DataFrame(iter=[1, 2], emax=[-0.1, -0.1])
        dfs = [empty_df(), good_df()]

        # ns_outputs / N_values length mismatch
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            dfs, [0, 1, 2], 4, μs, Ts)
        # N_values must include 0
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            dfs, [1, 2], 4, μs, Ts)
        # N_values must lie in 0:n_sites
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            dfs, [0, 5], 4, μs, Ts)
        # live_emax length mismatch
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            dfs, [0, 1], 4, μs, Ts; live_emax=[Float64[]])
        # a nonzero-N sector with no discarded samples and no live tail
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [empty_df(), empty_df()], [0, 1], 4, μs, Ts)
        # duplicate entries in N_values
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [empty_df(), good_df(), good_df()], [0, 1, 1], 4, μs, Ts)
        # non-positive temperature
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            dfs, [0, 1], 4, μs, [0.0u"K"])
    end

    # ================================================================
    @testset "_fixed_N_log_evidence: single-configuration sectors" begin
        # An empty ladder plus a live tail carries the sector's entire prior
        # mass, exactly 1 — independent of ω0 and of how many live energies
        # are supplied — so log_Z_NS = -βE exactly.
        fe = FreeBird.AnalysisTools._fixed_N_log_evidence
        empty_df = DataFrame(iter=Int[], emax=Float64[])
        Ts = [200.0, 400.0] .* u"K"
        E = -0.3
        K = 25
        for ω0 in (1.0, (K + 1) / K), live in (fill(E, K), [E])
            logZ, meanE = fe(empty_df, Ts;
                n_walkers=K, n_cull=1, ω0=ω0, live_energies=live, kb=kb)
            for (j, T) in enumerate(Ts)
                β = 1.0 / (kb * ustrip(u"K", T))
                @test isapprox(logZ[j], -β * E, rtol=1e-12)
                @test isapprox(meanE[j], E, rtol=1e-12)
            end
        end
    end

    # ================================================================
    @testset "Langmuir closed form from synthetic ladders" begin
        # On-site-only lattice gas: every N-particle configuration has
        # E = N ε_ads, so each per-N NS ladder is exactly flat and can be
        # synthesized without running a sampler. With Skilling weights
        # ω0 = (K+1)/K plus the live-set tail, Σω = 1 per sector (up to
        # O(r^n/K)) and the assembly must reproduce the Langmuir closed form
        # Ξ = (1 + z e^{-βε})^M to floating-point accuracy.
        M = 16
        ε = -0.04
        K = 25
        n_iters = 600
        ω0 = (K + 1) / K
        N_values = collect(0:M)

        dfs = [DataFrame(iter=collect(1:n_iters), emax=fill(N * ε, n_iters))
               for N in N_values]
        live = [fill(N * ε, K) for N in N_values]
        # The N = M sector uses the documented single-configuration recipe
        # (empty DataFrame + live tail), pinning it at tight tolerance: the
        # tail carries exactly the sector's prior mass 1, matching the flat
        # ladders' 1 + O(r^n/K) closure within rtol.
        dfs[end] = DataFrame(iter=Int[], emax=Float64[])

        μ_grid = [-0.08, -0.04, 0.0] .* u"eV"
        T_grid = [250.0, 300.0, 400.0] .* u"K"

        stats = gc_thermodynamic_stats_fixed_N(dfs, N_values, M, μ_grid, T_grid;
            n_walkers=K, n_cull=1, ω0=ω0, live_emax=live)

        @test stats.logXi isa Matrix{Float64}
        @test size(stats.logXi) == (length(μ_grid), length(T_grid))

        for (j, T) in enumerate(T_grid), (k, μ) in enumerate(μ_grid)
            β = 1.0 / (kb * ustrip(u"K", T))
            zeff = exp(β * (ustrip(u"eV", μ) - ε))
            θ = zeff / (1 + zeff)
            @test isapprox(stats.logXi[k, j], M * log1p(zeff), rtol=1e-8)
            @test isapprox(stats.mean_N[k, j], M * θ, rtol=1e-8)
            @test isapprox(stats.var_N[k, j], M * θ / (1 + zeff), rtol=1e-8)
            @test isapprox(stats.mean_U[k, j], ε * M * θ, rtol=1e-8)
        end
    end

    # ================================================================
    @testset "exact-enumeration validation (3x3 lattice, real NS runs)" begin
        using Random
        # The NS tolerance checks below are statistical: seed the global RNG
        # so they do not depend on which test files ran before this one.
        # (The random_seed field on the NS parameters is not consumed.)
        Random.seed!(1000)

        M = 9
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

        # A 3x3 periodic square lattice with the first N sites occupied.
        function lattice_at(N::Int)
            lat = MLattice{1,SquareLattice}(
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(3, 3, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1, 1.5],
                components=:equal,
                adsorptions=:full)
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            return lat
        end

        # Exact per-N energy lists from fixed-N enumeration (512 configs total).
        exact_E = Dict{Int,Vector{Float64}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(lattice_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            @test length(exact_E[N]) == binomial(M, N)
        end

        # Exact grand-canonical reference from the full spectrum.
        function exact_stats(μ_val::Float64, T_val::Float64)
            β = 1.0 / (kb * T_val)
            log_terms = Float64[]
            Ns = Int[]
            Es = Float64[]
            for N in 0:M, E in exact_E[N]
                push!(log_terms, N * β * μ_val - β * E)
                push!(Ns, N)
                push!(Es, E)
            end
            max_log = maximum(log_terms)
            w = exp.(log_terms .- max_log)
            sw = sum(w)
            logXi = max_log + log(sw)
            mean_N = sum(w .* Ns) / sw
            var_N = sum(w .* (Ns .^ 2)) / sw - mean_N^2
            mean_U = sum(w .* Es) / sw
            return logXi, mean_N, var_N, mean_U
        end

        # Canonical NS at each intermediate N; the walk (site hops) conserves N.
        K = 100
        n_iters = 800
        save_strategy = SaveEveryN(
            df_filename="test_lfn_df.csv",
            wk_filename="test_lfn.traj.extxyz",
            ls_filename="test_lfn.ls.extxyz",
            n_traj=1_000_000, n_snap=1_000_000, n_info=1_000_000)

        dfs = Vector{DataFrame}(undef, M + 1)
        live = Vector{Vector{Float64}}(undef, M + 1)

        # N = 0: empty lattice, handled internally (DataFrame ignored).
        dfs[1] = DataFrame(iter=Int[], emax=Float64[])
        live[1] = Float64[]

        for N in 1:(M - 1)
            walkers = [
                begin
                    lat = lattice_at(N)
                    generate_random_new_lattice_sample!(lat)
                    LatticeWalker(lat)
                end for _ in 1:K]
            liveset = LatticeGasWalkers(walkers, ham, perturb_energy=1e-12)
            ns_params = LatticeNestedSamplingParameters(
                mc_steps=60,
                allowed_fail_count=100_000)
            df, final_ls, _ = nested_sampling(
                liveset, ns_params, n_iters, MCRandomWalkClone(), save_strategy)
            dfs[N + 1] = df
            live[N + 1] = [ustrip(u"eV", w.energy) for w in final_ls.walkers]
        end

        # N = M: a single configuration, so NS cannot progress; the live-set
        # tail carries the whole sector (empty DataFrame + K copies of E_full).
        dfs[M + 1] = DataFrame(iter=Int[], emax=Float64[])
        live[M + 1] = fill(only(exact_E[M]), K)

        rm("test_lfn_df.csv", force=true)
        rm("test_lfn.traj.extxyz", force=true)
        rm("test_lfn.ls.extxyz", force=true)

        μ_grid = [-0.10, -0.06, -0.03] .* u"eV"
        T_grid = [300.0, 500.0] .* u"K"

        stats = gc_thermodynamic_stats_fixed_N(dfs, collect(0:M), M,
            μ_grid, T_grid;
            n_walkers=K, n_cull=1, ω0=(K + 1) / K, live_emax=live)

        for (j, T) in enumerate(T_grid), (k, μ) in enumerate(μ_grid)
            ref_logXi, ref_mean_N, ref_var_N, ref_mean_U = exact_stats(
                ustrip(u"eV", μ), ustrip(u"K", T))
            @test isapprox(stats.logXi[k, j], ref_logXi, atol=0.1)
            @test isapprox(stats.mean_N[k, j], ref_mean_N, rtol=0.05, atol=0.1)
            @test isapprox(stats.var_N[k, j], ref_var_N, rtol=0.25, atol=0.2)
            @test isapprox(stats.mean_U[k, j], ref_mean_U, rtol=0.05, atol=0.02)
        end

        # Monotonicity: ⟨N⟩ increases with μ at fixed T.
        for j in eachindex(T_grid)
            @test issorted(stats.mean_N[:, j])
        end
    end

end
