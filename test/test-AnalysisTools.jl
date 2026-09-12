@testset "AnalysisTools.jl tests" begin
    @testset "ωᵢ function tests" begin
        # Basic tests
        @test ωᵢ([1], 4) ≈ [1/5]
        @test ωᵢ([1, 2, 3], 4) ≈ [1/5, 1/5 * (4/5), 1/5 * (4/5)^2]
        @test ωᵢ([100], 4)[1] ≈ 1/5 * (4/5)^99 rtol=1e-10
        
        # Edge cases and properties
        @test ωᵢ(Int[], 4) == Float64[]
        @test_throws MethodError ωᵢ([1.5], 4)
        result = ωᵢ(collect(1:5), 4)
        @test all(diff(result) .< 0)  # Monotonic decrease
        @test all(result .> 0)        # All positive
    end

    @testset "log_ωᵢ function tests" begin
        # Agrees with log.(ωᵢ(...)) everywhere ωᵢ has not underflowed
        @test log_ωᵢ([1], 4) ≈ log.(ωᵢ([1], 4))
        @test log_ωᵢ(collect(1:200), 4) ≈ log.(ωᵢ(collect(1:200), 4))
        @test log_ωᵢ([1, 2, 3], 4; ω0=2.5) ≈ log.(ωᵢ([1, 2, 3], 4; ω0=2.5))
        @test log_ωᵢ([5], 7; n_cull=3) ≈ log.(ωᵢ([5], 7; n_cull=3))

        # ... and stays finite where it has. This is the whole point: ωᵢ is a
        # factor below one raised to the iteration number, so it reaches exactly
        # 0.0 and log() of it is -Inf, silently dropping the deepest samples
        # from any log-sum-exp.
        @test ωᵢ([100_000], 4)[1] == 0.0
        @test isinf(log(ωᵢ([100_000], 4)[1]))
        @test isfinite(log_ωᵢ([100_000], 4)[1])
        @test log_ωᵢ([100_000], 4)[1] ≈ log(1/5) + 99_999 * log(4/5)

        # Same contracts as ωᵢ
        @test log_ωᵢ(Int[], 4) == Float64[]
        @test_throws MethodError log_ωᵢ([1.5], 4)
        @test all(diff(log_ωᵢ(collect(1:5), 4)) .< 0)  # Monotonic decrease
    end

    @testset "gc_thermodynamic_stats live-set tail" begin
        # Nested sampling stops after finitely many iterations; the prior volume
        # left over, X_f = (K/(K+C))^{i_f}, is carried by the K survivors. These
        # tests pin the weight the tail gets, exactly.
        K, C = 4, 1
        r = K / (K + C)                      # 4/5
        n = 3
        iters = collect(1:n)

        # The default ω0 = 1 makes the dead-sample weights sum to 1 − X_f, so
        # with the tail the weights close to exactly 1. Tag dead samples with N = 0
        # and the live ones with N = 1: at β = 0 the returned ⟨N⟩ *is* the
        # fraction of the total weight carried by the live set, which must be
        # X_f = r^n. Nothing else in the pipeline can produce that number.
        ω0 = 1.0
        df = DataFrame(iter = iters,
                       omega = fill(1.0, n),
                       energy = fill(1.0, n),
                       num_particles = fill(0, n))
        live_E = fill(1.0, K)
        live_N = fill(1, K)

        _, _, meanN = gc_thermodynamic_stats(df, [0.0], K, 0.0;
                                             n_cull=C, ω0=ω0,
                                             live_energies=live_E,
                                             live_numbers=live_N)
        @test meanN[1] ≈ r^n

        # Without the live set the same call returns the truncated estimate,
        # which sees no particles at all.
        _, _, meanN_trunc = gc_thermodynamic_stats(df, [0.0], K, 0.0;
                                                   n_cull=C, ω0=ω0)
        @test meanN_trunc[1] ≈ 0.0

        # The bias the tail removes is low-temperature-first. Recorded samples
        # sit at E = 1, survivors at E = 0; as β grows the survivors must take
        # over, and do not if the tail is dropped.
        df2 = DataFrame(iter = iters,
                        omega = fill(1.0, n),
                        energy = fill(1.0, n),
                        num_particles = fill(0, n))
        meanE_tail, = gc_thermodynamic_stats(df2, [50.0], K, 0.0;
                                             n_cull=C, ω0=ω0,
                                             live_energies=zeros(K),
                                             live_numbers=zeros(Int, K))
        meanE_trunc, = gc_thermodynamic_stats(df2, [50.0], K, 0.0;
                                              n_cull=C, ω0=ω0)
        @test meanE_tail[1] ≈ 0.0 atol=1e-8
        @test meanE_trunc[1] ≈ 1.0

        # Ω for the live set is derived as E − μN, so a live set recorded
        # without a μ can still be reweighted at any μ.
        _, _, meanN_mu = gc_thermodynamic_stats(df, [0.0], K, -0.5;
                                                n_cull=C, ω0=ω0,
                                                live_energies=live_E,
                                                live_numbers=live_N)
        @test meanN_mu[1] ≈ r^n   # β = 0, so μ cannot matter here

        # Argument validation
        @test_throws ArgumentError gc_thermodynamic_stats(df, [0.0], K, 0.0;
                                                          live_energies=live_E)
        @test_throws DimensionMismatch gc_thermodynamic_stats(df, [0.0], K, 0.0;
                                                              live_energies=live_E,
                                                              live_numbers=[1])
    end

    @testset "gc_thermodynamic_stats evaluates Ω at the requested μ" begin
        # (Ω, E, N) are all recorded so a run can be reweighted to a μ it was
        # not sampled at. That only works if Ω is rebuilt from E and N at the
        # requested μ; reusing the recorded omega column would weight the
        # samples on the run's Ω-ladder while correcting Cv at the new μ.
        K, C = 4, 1
        iters = [1, 2]
        μ_run = -0.5
        Es = [0.0, -1.0]
        Ns = [0, 1]
        df = DataFrame(iter = iters,
                       omega = Es .- μ_run .* Ns,
                       energy = Es,
                       num_particles = Ns)
        β = 2.0
        w = exp.(log_ωᵢ(iters, K; n_cull=C))

        # At the run's own μ, nothing changes: derived Ω == recorded Ω.
        _, _, n_run = gc_thermodynamic_stats(df, [β], K, μ_run; n_cull=C)
        t_run = w .* exp.(-β .* (Es .- μ_run .* Ns))
        @test n_run[1] ≈ sum(t_run .* Ns) / sum(t_run)

        # At a different μ, Ω is rebuilt ...
        _, _, n_rw = gc_thermodynamic_stats(df, [β], K, 0.0; n_cull=C)
        t_rw = w .* exp.(-β .* Es)
        @test n_rw[1] ≈ sum(t_rw .* Ns) / sum(t_rw)

        # ... and that is a genuinely different answer from reusing df.omega,
        # so this test would fail if the column were read back.
        t_stale = w .* exp.(-β .* df.omega)
        @test !isapprox(n_rw[1], sum(t_stale .* Ns) / sum(t_stale); rtol=1e-3)
    end

    @testset "gc_thermodynamic_stats reports C_E, C_Ω and C_N" begin
        K, C = 4, 1
        iters = collect(1:6)
        Es = [0.0, -0.5, -1.0, -1.2, -1.5, -2.0]
        Ns = [0, 1, 1, 2, 2, 3]
        β = 3.0

        mkdf(μ) = DataFrame(iter = iters, omega = Es .- μ .* Ns,
                            energy = Es, num_particles = Ns)

        # The first three fields, in order, are still (mean_E, cv, mean_N), so
        # positional destructuring keeps working.
        r = gc_thermodynamic_stats(mkdf(-0.4), [β], K, -0.4; n_cull=C)
        a, b, c = gc_thermodynamic_stats(mkdf(-0.4), [β], K, -0.4; n_cull=C)
        @test a == r.mean_E
        @test b == r.cv
        @test c == r.mean_N
        @test all(isfinite, (r.cv[1], r.c_omega[1], r.c_N[1], r.var_N[1]))
        @test r.var_N[1] >= 0.0

        # C_Ω − C_E = −μ(∂⟨N⟩/∂T)_μ, so the two coincide at μ = 0 and only
        # there. This is the identity that makes reporting both meaningful
        # rather than redundant.
        r0 = gc_thermodynamic_stats(mkdf(0.0), [β], K, 0.0; n_cull=C)
        @test r0.c_omega[1] ≈ r0.cv[1]

        # And away from μ = 0, C_Ω is checked against a weighted Var(Ω) built
        # independently here — an exact check, rather than asserting only that
        # the two definitions happen to differ, which depends on the fixture.
        μt = -0.4
        kb = 8.617333262e-5
        rμ = gc_thermodynamic_stats(mkdf(μt), [β], K, μt; n_cull=C)
        Ω = Es .- μt .* Ns
        lt = log_ωᵢ(iters, K; n_cull=C) .- β .* Ω
        wts = exp.(lt .- maximum(lt))
        wts ./= sum(wts)
        varΩ = sum(wts .* Ω .^ 2) - sum(wts .* Ω)^2
        @test rμ.c_omega[1] ≈ kb * β^2 * varΩ

        # With N pinned, the particle-number projection is undefined rather than
        # zero, and all three definitions collapse onto k_Bβ²Var(E).
        df_fixedN = DataFrame(iter = iters, omega = Es .- (-0.4) .* fill(2, 6),
                              energy = Es, num_particles = fill(2, 6))
        rf = gc_thermodynamic_stats(df_fixedN, [β], K, -0.4; n_cull=C)
        @test rf.var_N[1] ≈ 0.0 atol=1e-12
        @test rf.c_N[1] ≈ rf.cv[1]
        @test rf.c_omega[1] ≈ rf.cv[1]

        # The empty result carries the same fields rather than a shorter tuple.
        rn = gc_thermodynamic_stats(1.0, Float64[], Float64[], Float64[], Int[], 0.0)
        @test all(isnan, (rn.mean_E, rn.cv, rn.mean_N, rn.c_omega, rn.c_N, rn.var_N))
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
end
