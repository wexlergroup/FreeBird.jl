# Effective-sample-size reductions for ideal-gas-referenced ledgers:
# gc_effective_sample_size_ideal_ref (modes :kish / :anchored /
# :anchored_uniform, relative normalization, explicit T0), welded to the
# N_eff field of gc_thermodynamic_stats_ideal_ref and to closed forms.
@testset "gc_effective_sample_size_ideal_ref" begin
    using Random

    kb_ess = 8.617333262e-5

    # Shared lattice fixture: seeded synthetic igref ledger with energies and
    # a live tail (descending emax, dense iter, K = 20)
    rng_ess = MersenneTwister(91201)
    n_fix = 400
    df_fix = DataFrame(iter=1:n_fix,
                       emax=sort(randn(rng_ess, n_fix) .* 0.05; rev=true),
                       num_particles=rand(rng_ess, 0:12, n_fix))
    live_E_fix = sort(randn(rng_ess, 20) .* 0.01; rev=true) .+
                 (minimum(df_fix.emax) - 0.05)
    live_N_fix = rand(rng_ess, 0:12, 20)
    K_fix = 20
    mus_fix = [-0.08, -0.02, 0.0, 0.02]
    Ts_fix = [200.0, 400.0]
    z0_fix = 0.7

    # Contrast fixture: athermal, N = 0 on the shallow half and N = 8 on the
    # deep half of a 4-nat ladder, so a positive reweighting exponent can
    # counteract the geometric shell decay inside the mu window
    n_con, K_con = 2000, 500
    df_con = DataFrame(iter=1:n_con, emax=zeros(n_con),
                       num_particles=vcat(zeros(Int, 1000), fill(8, 1000)))
    r_con = K_con / (K_con + 1)
    closed_con = (1 + r_con) * (1 - r_con^n_con)^2 /
                 ((1 - r_con) * (1 - r_con^(2 * n_con)))

    @testset ":kish welds to the N_eff field (lattice)" begin
        for (le, ln_) in ((nothing, nothing), (live_E_fix, live_N_fix))
            ess = gc_effective_sample_size_ideal_ref(
                df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix;
                live_emax=le, live_numbers=ln_)
            st = gc_thermodynamic_stats_ideal_ref(
                df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix;
                live_emax=le, live_numbers=ln_)
            @test size(ess) == (4, 2)
            @test all(isapprox.(ess, st.N_eff; rtol=1e-12))
        end
    end

    @testset "anchored closed form and run-independence at the anchor" begin
        # z0 = 1 makes the reference chemical potential exactly zero, so the
        # grid point mu = 0.0 evaluates the reweighting exponent as exact 0.0
        rng_b = MersenneTwister(91211)
        df_b = DataFrame(iter=1:n_con, emax=sort(rand(rng_b, n_con); rev=true),
                         num_particles=rand(rng_b, 0:16, n_con))
        eA = gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, [0.0], [300.0], K_con; mode=:anchored)
        eB = gc_effective_sample_size_ideal_ref(
            df_b, 16, 1.0, [0.0], [300.0], K_con; mode=:anchored)
        @test isapprox(eA[1, 1], closed_con; rtol=1e-12)
        # Two structurally different fixtures (different energies and particle
        # numbers, same J and K): bitwise-equal anchor values
        @test eA[1, 1] == eB[1, 1]
    end

    @testset ":anchored_uniform anchor, monotone rays, :kish non-monotonicity" begin
        mus_con = collect(range(-0.03, 0.03; length=61))
        i0 = findfirst(==(0.0), mus_con)
        eu = gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, mus_con, [300.0], K_con; mode=:anchored_uniform)
        @test eu[i0, 1] == Float64(n_con)
        @test all(diff(eu[i0:end, 1]) .<= 1e-12)
        @test all(diff(reverse(eu[1:i0, 1])) .<= 1e-12)
        # On the same fixture the raw Kish mode rises well above its anchor
        # value along the ray (measured peak ratio 1.572 at mu = 0.006): the
        # anchor is a baseline, not a maximum
        ek = gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, mus_con, [300.0], K_con)
        @test maximum(ek[:, 1]) > 1.2 * ek[i0, 1]
        # With all energies zero the default :anchored mode is bitwise the
        # same surface, so the non-monotonicity statement covers both Kish
        # modes on this fixture
        @test gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, mus_con, [300.0], K_con; mode=:anchored) == ek
    end

    @testset "explicit T0 anchoring" begin
        df_t = DataFrame(iter=1:5, emax=[0.4, 0.3, 0.2, 0.1, 0.05],
                         num_particles=[1, 2, 3, 2, 4])
        K_t, T_t, T0_t, mu_t, z0_t = 3, 350.0, 250.0, -0.03, 0.6
        beta_t = 1.0 / (kb_ess * T_t)
        beta0_t = 1.0 / (kb_ess * T0_t)
        s_t = beta_t * mu_t - log(z0_t)
        lw_hand = [log(1 / (K_t + 1)) + q * log(K_t / (K_t + 1)) +
                   df_t.num_particles[q] * s_t -
                   (beta_t - beta0_t) * df_t.emax[q] for q in 1:5]
        w_hand = exp.(lw_hand .- maximum(lw_hand))
        e_anc = gc_effective_sample_size_ideal_ref(
            df_t, 16, z0_t, [mu_t], [T_t], K_t; mode=:anchored, T0=T0_t)
        @test isapprox(e_anc[1, 1], sum(w_hand)^2 / sum(w_hand .^ 2); rtol=1e-12)
        lr_hand = [df_t.num_particles[q] * s_t -
                   (beta_t - beta0_t) * df_t.emax[q] for q in 1:5]
        r_hand = exp.(lr_hand .- maximum(lr_hand))
        e_unif = gc_effective_sample_size_ideal_ref(
            df_t, 16, z0_t, [mu_t], [T_t], K_t; mode=:anchored_uniform, T0=T0_t)
        @test isapprox(e_unif[1, 1], sum(r_hand)^2 / sum(r_hand .^ 2); rtol=1e-12)
        # T0 equal to the target temperature reproduces the default
        # per-target path bitwise (beta0 evaluates through the identical
        # expression, so the anchored exponent is exact 0.0 both ways)
        e_def = gc_effective_sample_size_ideal_ref(
            df_t, 16, z0_t, [mu_t], [T_t], K_t; mode=:anchored)
        e_T0T = gc_effective_sample_size_ideal_ref(
            df_t, 16, z0_t, [mu_t], [T_t], K_t; mode=:anchored, T0=T_t)
        @test e_T0T == e_def
    end

    @testset "athermal ledger: T-independence at mu = 0 only" begin
        rng_a = MersenneTwister(91221)
        n_a, K_a = 300, 30
        df_a = DataFrame(iter=1:n_a, emax=zeros(n_a),
                         num_particles=rand(rng_a, 0:10, n_a))
        ek = gc_effective_sample_size_ideal_ref(
            df_a, 16, 0.5, [0.0, 0.05], [200.0, 400.0], K_a)
        # At mu = 0 the exponent collapses to -ln z0 for every T: bitwise
        # T-independent. At mu != 0 the beta*mu term survives and the columns
        # legitimately differ.
        @test ek[1, 1] == ek[1, 2]
        @test ek[2, 1] != ek[2, 2]
        # With all energies exactly zero, :kish coincides bitwise with
        # :anchored under the default per-target T0
        ea = gc_effective_sample_size_ideal_ref(
            df_a, 16, 0.5, [0.0, 0.05], [200.0, 400.0], K_a; mode=:anchored)
        @test ek == ea
    end

    @testset "omega0 scale invariance (dead points only)" begin
        for mode in (:kish, :anchored, :anchored_uniform)
            e1 = gc_effective_sample_size_ideal_ref(
                df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix; mode=mode)
            e2 = gc_effective_sample_size_ideal_ref(
                df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix; mode=mode, ω0=1.05)
            e3 = gc_effective_sample_size_ideal_ref(
                df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix; mode=mode, ω0=21.0)
            @test all(isapprox.(e1, e2; rtol=1e-12))
            @test all(isapprox.(e1, e3; rtol=1e-12))
        end
    end

    @testset "relative normalization" begin
        mus_r = [-0.02, 0.0, 0.03]
        Ts_r = [250.0, 350.0]
        e_abs = gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, mus_r, Ts_r, K_con)
        e_rel = gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, mus_r, Ts_r, K_con; relative=true)
        @test all(isapprox.(e_rel[2, :], 1.0; rtol=1e-12))
        for j in 1:2
            @test all(isapprox.(e_rel[:, j], e_abs[:, j] ./ e_abs[2, j]; rtol=1e-12))
        end
        eur = gc_effective_sample_size_ideal_ref(
            df_con, 16, 1.0, mus_r, Ts_r, K_con; mode=:anchored_uniform, relative=true)
        @test all(eur[2, :] .== 1.0)
    end

    @testset "symmetric degradation on a particle-hole-symmetric fixture" begin
        # z0 = 1, E = 0, N ~ iid Binomial(16, 1/2): degradation is symmetric
        # in ln(z/z0) in distribution. Gate calibrated from seeds 91101-91103
        # at >= 3x the maximum |ratio - 1| = 0.0764, which the shipped seed
        # 91101 itself attains (gate 0.23 = 3.0x; seeds 91102/91103 measured
        # 0.0479/0.0283). Both flanks must also sit well below the anchored
        # baseline, so the symmetry claim is about a real degradation
        # (measured 0.637/0.689 of baseline at the shipped seed; gate 0.8).
        binom_ess(rng, M, p) = count(_ -> rand(rng) < p, 1:M)
        n_s, K_s = 3000, 1000
        rng_s = MersenneTwister(91101)
        df_s = DataFrame(iter=1:n_s, emax=zeros(n_s),
                         num_particles=[binom_ess(rng_s, 16, 0.5) for _ in 1:n_s])
        mu_s = 0.3 * kb_ess * 300.0
        e_s = gc_effective_sample_size_ideal_ref(
            df_s, 16, 1.0, [-mu_s, mu_s], [300.0], K_s)
        @test abs(e_s[2, 1] / e_s[1, 1] - 1) < 0.23
        r_s = K_s / (K_s + 1)
        base_s = (1 + r_s) * (1 - r_s^n_s)^2 / ((1 - r_s) * (1 - r_s^(2 * n_s)))
        @test e_s[1, 1] < 0.8 * base_s
        @test e_s[2, 1] < 0.8 * base_s
    end

    @testset "atomistic method" begin
        n_at, K_at = 150, 15
        rng_at = MersenneTwister(91231)
        df_at = DataFrame(iter=1:n_at,
                          emax=sort(rand(rng_at, n_at) .* 0.5; rev=true),
                          num_particles=rand(rng_at, 0:6, n_at),
                          log_compression=fill(log(K_at / (K_at + 1)), n_at))
        V_at = 1000.0u"Å^3"
        m_at = 39.948u"u"
        act_at = 0.01u"Å^-3"
        mu_at = [-0.15u"eV", -0.05u"eV"]
        T_at = [250.0u"K", 350.0u"K"]
        live_E_at = fill(0.0, 10)
        live_N_at = fill(3, 10)
        for (le, ln_) in ((nothing, nothing), (live_E_at, live_N_at))
            ess = gc_effective_sample_size_ideal_ref(
                df_at, V_at, m_at, act_at, mu_at, T_at;
                live_emax=le, live_numbers=ln_)
            st = gc_thermodynamic_stats_ideal_ref(
                df_at, V_at, m_at, act_at, mu_at, T_at;
                live_emax=le, live_numbers=ln_)
            @test all(isapprox.(ess, st.N_eff; rtol=1e-12))
        end
        # Zero-row live-only ledger is legal (mirrors the sibling's contract:
        # only ledgers with rows must carry :log_compression)
        df_empty_at = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[])
        e_live = gc_effective_sample_size_ideal_ref(
            df_empty_at, V_at, m_at, act_at, mu_at, T_at;
            live_emax=live_E_at, live_numbers=live_N_at)
        @test all(isfinite, e_live) && all(e_live .>= 1.0)
        # Unitful T0 dispatch on the anchored mode
        e_t0 = gc_effective_sample_size_ideal_ref(
            df_at, V_at, m_at, act_at, mu_at, T_at; mode=:anchored, T0=300.0u"K")
        @test size(e_t0) == (2, 2) && all(isfinite, e_t0)
    end

    @testset "argument validation" begin
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix; T0=300.0)
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix; mode=:bogus)
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix; live_emax=live_E_fix)
        @test_throws DimensionMismatch gc_effective_sample_size_ideal_ref(
            df_fix, 16, z0_fix, mus_fix, Ts_fix, K_fix;
            live_emax=live_E_fix, live_numbers=[1])
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            DataFrame(iter=Int[], emax=Float64[], num_particles=Int[]),
            16, z0_fix, mus_fix, Ts_fix, K_fix)
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            df_fix, 16, -1.0, mus_fix, Ts_fix, K_fix)
        # Atomistic mirrors: rows without :log_compression; corrupted column
        df_nolc = DataFrame(iter=1:3, emax=[0.3, 0.2, 0.1], num_particles=[1, 1, 1])
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            df_nolc, 1000.0u"Å^3", 39.948u"u", 0.01u"Å^-3",
            [-0.1u"eV"], [300.0u"K"])
        df_badlc = DataFrame(iter=1:3, emax=[0.3, 0.2, 0.1],
                             num_particles=[1, 1, 1],
                             log_compression=[-0.1, 0.5, -0.1])
        @test_throws ArgumentError gc_effective_sample_size_ideal_ref(
            df_badlc, 1000.0u"Å^3", 39.948u"u", 0.01u"Å^-3",
            [-0.1u"eV"], [300.0u"K"])
    end
end
