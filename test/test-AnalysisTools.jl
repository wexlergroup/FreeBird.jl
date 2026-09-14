@testset "AnalysisTools.jl tests" begin
    @testset "ωᵢ function tests" begin
        # Basic tests
        @test ωᵢ([1], 4) ≈ [1/5 * (4/5)^1]
        @test ωᵢ([1, 2, 3], 4) ≈ [1/5 * (4/5)^1, 1/5 * (4/5)^2, 1/5 * (4/5)^3]
        @test ωᵢ([0], 4) ≈ [1/5]
        @test ωᵢ([100], 4)[1] ≈ 1/5 * (4/5)^100 rtol=1e-10
        
        # Edge cases and properties
        @test ωᵢ(Int[], 4) == Float64[]
        @test_throws MethodError ωᵢ([1.5], 4)
        result = ωᵢ(collect(1:5), 4)
        @test all(diff(result) .< 0)  # Monotonic decrease
        @test all(result .> 0)        # All positive
    end

    @testset "Partition function and internal energy tests" begin
        # Basic functionality
        @test partition_function(1.0, [1.0], [0.0]) ≈ 1.0
        @test internal_energy(1.0, [1.0], [0.0]) == 0.0

        # Temperature limits
        ωi = collect(0.40:-0.1:0.10)
        Ei = collect(0.0:0.3:1.0)
        @test partition_function(1e-10, ωi, Ei) ≈ 1.0 rtol=1e-8
        @test internal_energy(1e-10, ωi, Ei) ≈ 0.3 rtol=1e-8

        # Edge cases
        @test partition_function(1.0, Float64[], Float64[]) == 0.0
        @test isnan(internal_energy(1.0, Float64[], Float64[]))
        @test_throws DimensionMismatch partition_function(1.0, [0.5, 0.5], Ei)
    end

    @testset "Constant-volume heat capacity tests" begin
        kb = 8.617333262e-5
        ωi = [0.4, 0.3, 0.2, 0.1]
        Ei = [0.0, 0.3, 0.6, 0.9]
        dof = 3

        # Basic functionality and limits
        @test cv(1.0, [1.0], [0.0], 3) ≈ 3 * kb / 2.0 rtol=1e-8
        @test cv(1e-10, ωi, Ei, dof) ≈ dof * kb / 2.0 rtol=1e-8
        
        # DataFrame interface
        df = DataFrame(iter = [1, 2, 3], emax = [1.0, 1.5, 2.0])
        β = 1.0
        βs = [0.1, 1.0, 10.0]
        n_walkers = 4
        result = cv(df, βs, dof, n_walkers)
        @test length(result) == length(βs)
        @test all(.!isnan.(result))

        # Edge cases
        @test cv(1.0, ωi, Ei, 0) ≥ 0
        @test isnan(cv(β, Float64[], Float64[], dof))
        @test_throws DimensionMismatch cv(1.0, [0.5, 0.5], Ei, dof)
    end

    # Test cv function (Wang-Landau sampling)
    @testset "cv_wang_landau" begin
        Ts = [300.0, 400.0]
        dof = 3
        energy_bins = [1.0, 2.0, 3.0]
        entropy = [0.1, 0.2, 0.3]
        cvs = cv(Ts, dof, energy_bins, entropy)
        @test length(cvs) == 2
    end

    @testset "read dataframe" begin
        # Create a sample DataFrame
        df = DataFrame(iter = [1, 2, 3], emax = [1.0, 1.5, 2.0])
        # Write to CSV
        write_df("test_output_df.csv", df)
        # Write to Arrow
        write_df("test_output_df.arrow", df)
        # Read from CSV
        df_read_csv = read_output("test_output_df.csv")
        # Read from Arrow
        df_read_arrow = read_output("test_output_df.arrow")
        @test df == df_read_csv
        @test df == df_read_arrow
        # Clean up
        rm("test_output_df.csv")
        rm("test_output_df.arrow")
    end

    @testset "microcanonical inflection-point analysis" begin
        # Synthesize an NS ladder whose density of states has a prescribed caloric
        # temperature β(E)=dS/dE: build g=exp(∫β), its normalized CDF G, and invert
        # G at the NS prior volumes X_i=(K/(K+n_cull))^i so emax(i)=G⁻¹(X_i).
        # Also returns the (Eg, G) grid so tests can check S_vol = ln G exactly.
        function synth_ladder(βfun; Emin=0.0, Emax=10.0, K=200, n=4000, ng=20001, n_cull=1)
            Eg = collect(range(Emin, Emax; length=ng)); dE = Eg[2] - Eg[1]
            βg = βfun.(Eg)
            S = zeros(ng); for j in 2:ng; S[j] = S[j-1] + 0.5*(βg[j]+βg[j-1])*dE; end
            g = exp.(S .- maximum(S))
            G = zeros(ng); for j in 2:ng; G[j] = G[j-1] + 0.5*(g[j]+g[j-1])*dE; end
            G ./= G[end]
            c = K/(K+n_cull); Es = Float64[]
            for i in 1:n
                Xi = c^i
                if Xi <= G[1]; push!(Es, Eg[1]); continue; end
                if Xi >= G[end]; push!(Es, Eg[end]); continue; end
                jj = searchsortedfirst(G, Xi)
                push!(Es, Eg[jj-1] + (Xi-G[jj-1])/(G[jj]-G[jj-1])*(Eg[jj]-Eg[jj-1]))
            end
            return DataFrame(iter = collect(1:n), emax = Es), K, Eg, G
        end
        # linear interpolation of the CDF grid at energy x
        function cdf_at(Eg, G, x)
            jj = clamp(searchsortedfirst(Eg, x), 2, length(Eg))
            G[jj-1] + (x - Eg[jj-1])/(Eg[jj] - Eg[jj-1])*(G[jj] - G[jj-1])
        end
        Etr = 5.0; w = 0.8

        @testset "entropy/β recovery, no spurious transition" begin
            βfun(E) = 2.0 - 0.10*E - 0.005*E^2          # β>0, γ strictly monotonic
            df, K, Eg, G = synth_ladder(βfun)
            E, S = microcanonical_entropy(df, K)
            @test issorted(E)
            @test all(isfinite, S)
            # volume entropy is exact by construction: S_vol = ln X_i = ln G(E_i)
            Ev, Sv = microcanonical_entropy(df, K; kind=:volume)
            @test all(abs(Sv[k] - log(cdf_at(Eg, G, Ev[k]))) < 1e-6
                      for k in eachindex(Ev) if 0.5 < Ev[k] < 9.5)
            d = caloric_derivatives(df, K; max_order=2)
            mid = findall(e -> 3.0 < e < 7.0, d.E)
            @test !isempty(mid)
            @test all(abs(d.β[i] - βfun(d.E[i])) < 0.15 for i in mid)   # β recovered
            # no false positives in the resolved bulk (edges/low-T flattening excluded)
            @test isempty(inflection_transitions(df, K; kb=1.0, energy_window=(2.0, 9.0)))
        end

        @testset "n_cull ≠ 1 prior-volume factor" begin
            βfun(E) = 2.0 - 0.12*E
            df, K = synth_ladder(βfun; n_cull=4)
            d = caloric_derivatives(df, K; n_cull=4)
            mid = findall(e -> 3.0 < e < 7.0, d.E)
            @test !isempty(mid)
            @test all(abs(d.β[i] - βfun(d.E[i])) < 0.15 for i in mid)   # cc = ln(K/(K+4))
        end

        @testset "second-order transition (γ negative peak)" begin
            βfun(E) = 2.0 - 0.15*E + 0.10*tanh((E-Etr)/w)   # C/w=0.125 < 0.15 (concave)
            df, K = synth_ladder(βfun)
            ts = inflection_transitions(df, K; kb=1.0, max_order=2)
            @test length(ts) == 1                       # exactly one transition reported
            t = ts[1]
            @test abs(t.E_tr - Etr) < 1.0
            @test t.order == 2
            @test t.kind == :independent
            @test isapprox(t.T_tr, 1/βfun(Etr); rtol=0.25)              # T_tr ≈ 1/β(Etr)
        end

        @testset "first-order transition (β backbending)" begin
            βfun(E) = 2.0 - 0.15*E + 0.40*tanh((E-Etr)/w)   # C/w=0.50 > 0.15 → backbend
            df, K = synth_ladder(βfun)
            ts = inflection_transitions(df, K; kb=1.0, max_order=2)
            @test length(ts) == 1                       # exactly one transition reported
            t = ts[1]
            @test abs(t.E_tr - Etr) < 1.5
            @test t.order == 1
            @test t.kind == :independent
        end

        @testset "transition_convergence" begin
            βfun(E) = 2.0 - 0.15*E + 0.10*tanh((E-Etr)/w)
            df1, _ = synth_ladder(βfun; K=150)
            df2, _ = synth_ladder(βfun; K=300)
            res = transition_convergence([df1, df2], [150, 300]; kb=1.0)
            @test res isa Vector
            @test all(r -> haskey(r, :converged) && haskey(r, :T_by_K) && length(r.T_by_K) == 2, res)
            @test any(r -> r.converged, res)        # deterministic DOS → stable across K
            # negative paths: zero tolerance flags any finite drift as unconverged
            res0 = transition_convergence([df1, df2], [150, 300]; kb=1.0, tol=0.0)
            @test !isempty(res0) && all(r -> !r.converged, res0)
            # an impossibly tight matching window leaves the lower-K ladder unmatched
            resm = transition_convergence([df1, df2], [150, 300]; kb=1.0, match_tol=1e-12)
            @test !isempty(resm)
            @test all(r -> r.T_drift == Inf && !r.converged, resm)
            @test all(r -> ismissing(r.T_by_K[1]) && !ismissing(r.T_by_K[2]), resm)
        end

        @testset "argument validation" begin
            df, K = synth_ladder(E -> 2.0 - 0.12*E)
            @test_throws ArgumentError microcanonical_entropy(df, K; kind=:bogus)
            @test_throws ArgumentError caloric_derivatives(df, K; max_order=5)
            @test_throws ArgumentError inflection_transitions(df, K; max_order=0)
            @test_throws ArgumentError inflection_transitions(df, K; ground_trim=0.6)
            @test_throws ArgumentError transition_convergence([df], [K])
            @test_throws DimensionMismatch transition_convergence([df, df], [K])
            @test_throws ArgumentError microcanonical_entropy(DataFrame(iter=Int[], emax=Float64[]), K)
            @test_throws ArgumentError microcanonical_entropy(df, 0)
            @test_throws ArgumentError microcanonical_entropy(DataFrame(a=collect(1:100)), K)
            @test_throws ArgumentError microcanonical_entropy(df, K; halfwidth=1)
            @test_throws ArgumentError microcanonical_entropy(df, K; halfwidth=-3)
            @test_throws ArgumentError inflection_transitions(df, K; edge=0)
        end
    end

    @testset "kish_effective_sample_size" begin
        using Random

        # Exact closures: uniform log-weights return exactly n (the ratio
        # n^2/n is exact); the two-point value against an independent
        # linear-space hand computation
        for n in (1, 2, 7, 100)
            @test kish_effective_sample_size(zeros(n)) == Float64(n)
            @test kish_effective_sample_size(fill(3.7, n)) == Float64(n)
        end
        let a = 0.3, b = -1.2
            wa, wb = exp(a), exp(b)
            @test isapprox(kish_effective_sample_size([a, b]),
                           (wa + wb)^2 / (wa^2 + wb^2); rtol=1e-12)
        end

        # Shift invariance: generic shifts at rtol (the addition rounds each
        # entry before the max-shift can remove it, so bitwise equality there
        # would be rounding luck); bitwise on an exact-dyadic fixture whose
        # entries and shift are multiples of 2^-16 with bounded magnitude, so
        # every addition is exact and the max-shift restores the residuals
        rng = MersenneTwister(91021)
        lw = randn(rng, 200)
        base = kish_effective_sample_size(lw)
        for c in (pi, -17.3)
            @test isapprox(kish_effective_sample_size(lw .+ c), base; rtol=1e-12)
        end
        # An extreme shift costs precision in the inputs themselves (ulp(1e6)
        # is ~1.2e-10, rounding every entry before the max-shift can act), so
        # its gate reflects input rounding, not the ratio: measured 3.5e-12
        @test isapprox(kish_effective_sample_size(lw .+ 1e6), base; rtol=1e-9)
        lwd = Float64.(rand(rng, -2^18:2^18, 300)) ./ 2^16
        cshift = 7.0 * 2.0^-5
        @test kish_effective_sample_size(lwd .+ cshift) ==
              kish_effective_sample_size(lwd)

        # Bounds, with both equality cases hit exactly (uniform above; a
        # dominant weight whose competitors underflow to exact zeros here)
        for seed in (91031, 91032, 91033)
            lwr = 3 .* randn(MersenneTwister(seed), 500)
            e = kish_effective_sample_size(lwr)
            @test 1.0 <= e <= 500.0
        end
        @test kish_effective_sample_size([0.0, -800.0, -800.0]) == 1.0

        # Geometric shell weights (the nested-sampling prior profile): the
        # closed form (1 + r)(1 - r^J)^2 / ((1 - r)(1 - r^(2J))) with
        # r = K/(K + 1), and its deep-ladder limit 2K + 1
        for (K, J) in ((10, 200), (500, 50_000))
            r = K / (K + 1)
            lw_geo = (1:J) .* log(r) .- log(K + 1)
            closed = (1 + r) * (1 - r^J)^2 / ((1 - r) * (1 - r^(2J)))
            @test isapprox(kish_effective_sample_size(lw_geo), closed; rtol=1e-12)
        end
        let K = 500, J = 50_000
            r = K / (K + 1)
            lw_geo = (1:J) .* log(r) .- log(K + 1)
            @test isapprox(kish_effective_sample_size(lw_geo), 2K + 1; rtol=1e-10)
        end

        # Independent-draw laws for iid weights a^N: ESS/n concentrates on
        # exp(-λ(a-1)^2) for N ~ Poisson(λ) and on
        # ((1 + p(a-1))^2 / (1 + p(a^2-1)))^M for N ~ Binomial(M, p).
        # Gates calibrated from three seeds each (91001-91003, 91011-91013)
        # at >= 3x the maximum relative deviation, per statistic; shipped
        # seeds 91001/91011 (measured deviations 0.00006/0.00141 and
        # 0.00010/0.02173 against gates 0.0015/0.07 and 0.006/0.24).
        binom_draw(rng, M, p) = count(_ -> rand(rng) < p, 1:M)
        function pois_draw(rng, lam)
            L = exp(-lam)
            k = 0
            acc = rand(rng)
            while acc > L
                k += 1
                acc *= rand(rng)
            end
            return k
        end
        n_iid = 40_000
        let rng_b = MersenneTwister(91001)
            Ns = [binom_draw(rng_b, 64, 0.5) for _ in 1:n_iid]
            for (s, gate) in ((0.05, 0.0015), (0.2, 0.07))
                a = exp(s)
                pred = ((1 + 0.5 * (a - 1))^2 / (1 + 0.5 * (a^2 - 1)))^64
                meas = kish_effective_sample_size(s .* Ns) / n_iid
                @test abs(meas / pred - 1) < gate
            end
        end
        let rng_p = MersenneTwister(91011)
            Ns = [pois_draw(rng_p, 8.0) for _ in 1:n_iid]
            for (s, gate) in ((0.1, 0.006), (0.3, 0.24))
                a = exp(s)
                pred = exp(-8.0 * (a - 1)^2)
                meas = kish_effective_sample_size(s .* Ns) / n_iid
                @test abs(meas / pred - 1) < gate
            end
        end

        # Errors and zero weights: -Inf entries are legal zero weights that
        # drop out; NaN, +Inf, emptiness, and all-(-Inf) throw
        @test kish_effective_sample_size([-Inf, 0.0, 0.0]) == 2.0
        @test_throws ArgumentError kish_effective_sample_size(Float64[])
        @test_throws ArgumentError kish_effective_sample_size([1.0, NaN])
        @test_throws ArgumentError kish_effective_sample_size([1.0, Inf])
        @test_throws ArgumentError kish_effective_sample_size([-Inf, -Inf])
    end
end