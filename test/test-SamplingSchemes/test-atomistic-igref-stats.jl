# Atomistic ledger assembly for ideal-gas-referenced grand-canonical runs.
#
# Calibration ledger:
# - Cross-route agreement (igref reduction vs fixed-N stitching, dilute LJ gas,
#   K = 24, per-T μ grids solved from target zV in {4, 6}; igref seeds
#   {1,2,3,62621} against fixed-N stitching seeds {1001,1002,1003,63621}):
#   max devs logXi 0.341, mean_N 0.636, N_eff floor observed 17.3; gates ship at
#   >= 3x (1.05 and 2.0) with an N_eff floor assertion at 10.
# - The exact-closure testset needs no calibration: the compression ledger is
#   CONSTRUCTED to encode the truncated Poisson prior exactly, so every
#   grand-canonical identity is a floating-point weld, gated at 1e-11 relative.
@testset "atomistic ideal-gas-referenced ledger assembly tests" begin
    using Random

    kb = 8.617333262e-5
    Vq = 1000.0u"Å^3"
    V = 1000.0
    mass_ar = 39.948u"u"
    lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass_ar, T * u"K"))
    mu_for(zV, T) = (kb * T * (log(zV / V) + 3 * log(lam(T)))) * u"eV"

    logfact(n) = n == 0 ? 0.0 : sum(log, 1:n)
    poisson_pmf(lam0, n) = exp(n * log(lam0) - lam0 - logfact(n))

    # A compression ledger that ENCODES the truncated Poisson(z0V) prior exactly:
    # shell j carries mass p(N = j - 1) at energy zero, the tail walker carries
    # the (negligible) truncation remainder
    function poisson_ledger(z0V; nmax=60)
        p = [poisson_pmf(z0V, n) for n in 0:nmax]
        remainder = max(1.0 - sum(p), 1.0e-300)
        masses = vcat(p, remainder)
        # X[j] = prior mass from shell j onward, positive by construction
        # (a forward 1 - cumsum cancels to negative dust at deep truncation)
        X = reverse(cumsum(reverse(masses)))
        lc = [log(X[j+1] / X[j]) for j in 1:nmax+1]
        df = DataFrame(iter=collect(1:nmax+1), emax=zeros(nmax + 1),
                       num_particles=collect(0:nmax), log_compression=lc)
        return df, [0.0], [nmax + 1]     # live tail: the truncation remainder
    end

    @testset "exact closure: the assembly reproduces the grand-canonical closed forms" begin
        z0V = 3.0
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = poisson_ledger(z0V)
        Ts = [250.0, 300.0, 400.0]
        for T in Ts
            mus = [mu_for(1.0, T), mu_for(3.0, T), mu_for(7.0, T)]
            stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                                     live_emax=live_e, live_numbers=live_n)
            for (i, zV) in enumerate([1.0, 3.0, 7.0])
                @test abs(stats.logXi[i, 1] - zV) < 1e-11 * max(zV, 1.0)
                @test abs(stats.mean_N[i, 1] - zV) < 1e-11 * zV
                @test abs(stats.var_N[i, 1] - zV) < 1e-10 * zV
                @test abs(stats.mean_U[i, 1]) < 1e-14
            end
        end
    end

    @testset "shell weld against the exported ωᵢ decode" begin
        z0V = 3.0
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = poisson_ledger(z0V; nmax=20)
        T = 300.0
        mus = [mu_for(2.0, T)]
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                                 live_emax=live_e, live_numbers=live_n)
        # independent hand assembly through the exported linear-space decode
        w = ωᵢ(Vector{Float64}(df.log_compression); ω0=1.0)
        tail = exp(sum(df.log_compression)) / length(live_e)
        s = 2.0 / z0V
        hand = z0V + log(sum(w .* s .^ df.num_particles) + tail * s^live_n[1])
        @test abs(stats.logXi[1, 1] - hand) < 1e-12
    end

    @testset "reweighting window collapse is visible in N_eff" begin
        z0V = 3.0
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = poisson_ledger(z0V)
        T = 300.0
        # zV = 200 sits far beyond the ledger's N <= 60 support, so the reweighted
        # mass concentrates on the last rows; within the support (zV = 3) the exact
        # ledger reweights with a healthy effective sample size
        mus = [mu_for(3.0, T), mu_for(200.0, T)]
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                                 live_emax=live_e, live_numbers=live_n)
        @test stats.N_eff[1, 1] > 5.0
        @test stats.N_eff[2, 1] < 2.5
        @test stats.N_eff[1, 1] / stats.N_eff[2, 1] > 2.0
    end

    @testset "ledger-shape hardening" begin
        z0 = (3.0 / V)u"Å^-3"
        T = [300.0]u"K"
        mus = [mu_for(3.0, 300.0)]
        # ended-mid-eviction ledger: tail splits over the SURVIVING count
        lc = [log(8 / 9), log(8 / 9), log(8 / 9), log(7 / 8), log(6 / 7)]
        df = DataFrame(iter=1:5, emax=zeros(5), num_particles=[3, 2, 4, 3, 2],
                       log_compression=lc)
        live_e = zeros(6)
        live_n = [2, 3, 1, 4, 2, 3]
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, T;
                                                 live_emax=live_e, live_numbers=live_n)
        w = ωᵢ(lc; ω0=1.0)
        tail_each = exp(sum(lc)) / 6
        s = 3.0 / 3.0
        hand = 3.0 + log(sum(w .* s .^ df.num_particles) + tail_each * sum(s .^ live_n))
        @test abs(stats.logXi[1, 1] - hand) < 1e-12
        # total prior mass closes: shells + tail account for exactly 1
        @test abs(sum(w) + tail_each * 6 - 1.0) < 1e-14
        # a ledger with rows but no compression column is rejected
        df_bad = DataFrame(iter=1:2, emax=zeros(2), num_particles=[1, 2])
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_bad, Vq, mass_ar, z0, mus, T; live_emax=live_e, live_numbers=live_n)
        # a corrupted compression column (positive entry) is rejected
        df_pos = DataFrame(iter=1:2, emax=zeros(2), num_particles=[1, 2],
                           log_compression=[log(0.9), 0.1])
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_pos, Vq, mass_ar, z0, mus, T; live_emax=live_e, live_numbers=live_n)
        # zero-accept run: empty ledger with a live-set-only reduction is exact
        df_empty = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[],
                             log_compression=Float64[])
        stats_live = gc_thermodynamic_stats_ideal_ref(df_empty, Vq, mass_ar, z0, mus, T;
                                                      live_emax=zeros(4),
                                                      live_numbers=[2, 3, 3, 4])
        hand_live = 3.0 + log(sum((1.0) .^ [2, 3, 3, 4]) / 4)
        @test abs(stats_live.logXi[1, 1] - hand_live) < 1e-12
        # empty ledger with no live walkers is rejected
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_empty, Vq, mass_ar, z0, mus, T)
        # mismatched live vectors are rejected
        @test_throws DimensionMismatch gc_thermodynamic_stats_ideal_ref(
            df_empty, Vq, mass_ar, z0, mus, T; live_emax=zeros(3), live_numbers=[1, 2])
    end

    @testset "cross-route agreement on a dilute interacting gas (seed 62621)" begin
        box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        pbc = (true, true, true)
        seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
        mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                               empty(position(seed_at, :)), empty(species(seed_at, :)),
                               empty(mass(seed_at, :)))
        V12 = 1728.0
        Vq12 = 1728.0u"Å^3"
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
        z0 = (6.0 / V12)u"Å^-3"
        mu12(zV, T) = (kb * T * (log(zV / V12) + 3 * log(lam(T)))) * u"eV"
        Ts = [300.0, 400.0]
        mugrids = [[mu12(4.0, T), mu12(6.0, T)] for T in Ts]
        mksave(tag) = SaveEveryN(df_filename="_igstats_$(tag).csv",
                                 wk_filename="_igstats_$(tag).traj.extxyz",
                                 ls_filename="_igstats_$(tag).ls.extxyz",
                                 n_traj=10^7, n_snap=10^7, n_info=10^7)
        clean(tag) = for f in ["_igstats_$(tag).csv", "_igstats_$(tag).traj.extxyz", "_igstats_$(tag).ls.extxyz"]
            rm(f, force=true)
        end

        Random.seed!(62621)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:24], lj)
        params = AtomisticIGRefGCNSParameters(mc_steps=60, reference_activity=z0,
                                              species=:Ar, allowed_fail_count=100)
        df, lso, _ = ideal_gas_referenced_nested_sampling(
            ls, params, 300, MCAtomGrandCanonicalMoves(), mksave("ig"))
        clean("ig")
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]

        Random.seed!(63621)
        Nmax = 12
        dfs = DataFrame[]
        lives = Vector{Float64}[]
        for N in 0:Nmax
            if N == 0
                push!(dfs, DataFrame(iter=Int[], emax=Float64[]))
                push!(lives, zeros(8))
                continue
            end
            walkers = AtomWalker{1}[]
            for _ in 1:24
                coords = [[rand() * 12.0, rand() * 12.0, rand() * 12.0] for _ in 1:N]
                push!(walkers, AtomWalker{1}(FastSystem(atomic_system(
                    [:Ar => SVector{3}(c)u"Å" for c in coords], box, pbc))))
            end
            lsN = GenericAtomWalkers(walkers, lj)
            p = NestedSamplingParameters(mc_steps=60, initial_step_size=0.5, step_size=0.5,
                                         step_size_lo=0.01, step_size_up=2.0,
                                         allowed_fail_count=1000)
            dfN, lsoN, _ = nested_sampling(lsN, p, 150, MCRandomWalkClone(), mksave("fx$N"))
            clean("fx$N")
            push!(dfs, dfN)
            push!(lives, [ustrip(u"eV", w.energy) for w in lsoN.walkers])
        end

        for k in 1:2
            ig = gc_thermodynamic_stats_ideal_ref(df, Vq12, mass_ar, z0, mugrids[k], [Ts[k]]u"K";
                                                  live_emax=live_e, live_numbers=live_n)
            fx = gc_thermodynamic_stats_fixed_N(dfs, collect(0:Nmax), Vq12, mass_ar,
                                                mugrids[k], [Ts[k]]u"K";
                                                n_walkers=24, live_emax=lives)
            @test maximum(abs.(ig.logXi .- fx.logXi)) < 1.05
            @test maximum(abs.(ig.mean_N .- fx.mean_N)) < 2.0
            @test minimum(ig.N_eff) > 10.0
        end
    end
end

@testset "particle-number distributions from the reductions" begin
    kb = 8.617333262e-5
    Vq = 1000.0u"Å^3"
    V = 1000.0
    mass_ar = 39.948u"u"
    lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass_ar, T * u"K"))
    mu_for(zV, T) = (kb * T * (log(zV / V) + 3 * log(lam(T)))) * u"eV"
    logfact(n) = n == 0 ? 0.0 : sum(log, 1:n)
    poisson_pmf(lam0, n) = exp(n * log(lam0) - lam0 - logfact(n))
    function poisson_ledger(z0V; nmax=60)
        p = [poisson_pmf(z0V, n) for n in 0:nmax]
        remainder = max(1.0 - sum(p), 1.0e-300)
        masses = vcat(p, remainder)
        X = reverse(cumsum(reverse(masses)))
        lc = [log(X[j+1] / X[j]) for j in 1:nmax+1]
        df = DataFrame(iter=collect(1:nmax+1), emax=zeros(nmax + 1),
                       num_particles=collect(0:nmax), log_compression=lc)
        return df, [0.0], [nmax + 1]
    end

    @testset "exact Poisson closure of p_N" begin
        z0V = 3.0
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = poisson_ledger(z0V)
        T = 300.0
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0,
            [mu_for(2.0, T), mu_for(5.0, T)], [T]u"K";
            live_emax=live_e, live_numbers=live_n)
        @test stats.N_support == 0:61
        @test size(stats.p_N) == (2, 1, 62)
        for (i, zV) in enumerate([2.0, 5.0])
            for n in 0:15
                @test abs(stats.p_N[i, 1, n + 1] - poisson_pmf(zV, n)) < 1e-11
            end
            @test abs(sum(stats.p_N[i, 1, :]) - 1.0) < 1e-12
        end
    end

    @testset "moment identities per grid point" begin
        z0V = 4.0
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = poisson_ledger(z0V; nmax=40)
        Ts = [250.0, 350.0]
        mus = [mu_for(1.5, 300.0), mu_for(4.0, 300.0)]
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, Ts * u"K";
            live_emax=live_e, live_numbers=live_n)
        Ns = collect(stats.N_support)
        for i in 1:2, j in 1:2
            m1 = sum(stats.p_N[i, j, :] .* Ns)
            m2 = sum(stats.p_N[i, j, :] .* Ns .^ 2)
            @test isapprox(m1, stats.mean_N[i, j]; rtol=1e-12)
            # The two routes reassociate m2 - m1^2 differently, so the identity
            # carries an eps*m2 cancellation floor where the variance collapses
            @test isapprox(m2 - m1^2, stats.var_N[i, j]; rtol=1e-11,
                           atol=1e-11 * max(1.0, m2))
        end
    end

    @testset "live-tail mixture weld, bit-for-bit" begin
        z0V = 2.5
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = poisson_ledger(z0V; nmax=8)
        T = 300.0
        mu = mu_for(3.0, T)
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, [mu], [T]u"K";
            live_emax=live_e, live_numbers=live_n)
        # hand assembly replicating the source arithmetic exactly
        lc = Vector{Float64}(df.log_compression)
        cs = cumsum(lc)
        log_w0 = vcat(0.0 .+ vcat(0.0, cs[1:end-1]) .+ log.(-expm1.(lc)),
                      [cs[end] - log(1)])
        Es = vcat(zeros(9), live_e)
        Nsamp = vcat(Vector{Float64}(df.num_particles), Float64.(live_n))
        β = 1.0 / (kb * T)
        s = β * ustrip(u"eV", mu) - 3.0 * log(lam(T)) - log(ustrip(u"Å^-3", z0))
        log_w = log_w0 .+ s .* Nsamp .- β .* Es
        max_log = maximum(log_w)
        ws = exp.(log_w .- max_log)
        sum_w = sum(ws)
        acc = zeros(Float64, maximum(Int.(Nsamp)) + 1)
        for k in eachindex(ws)
            acc[Int(Nsamp[k]) + 1] += ws[k]
        end
        @test stats.p_N[1, 1, :] == acc ./ sum_w
    end

    @testset "support holes stay exactly zero" begin
        # masses on N in {0, 2} only, live tail at N = 4: N = 1 and N = 3 are
        # in-support holes that must carry exactly 0.0
        masses = [0.4, 0.35, 0.25]
        X = reverse(cumsum(reverse(masses)))
        lc = [log(X[j+1] / X[j]) for j in 1:2]
        df = DataFrame(iter=[1, 2], emax=zeros(2), num_particles=[0, 2],
                       log_compression=lc)
        T = 300.0
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, (2.0 / V)u"Å^-3",
            [mu_for(2.0, T)], [T]u"K"; live_emax=[0.0], live_numbers=[4])
        @test stats.N_support == 0:4
        @test stats.p_N[1, 1, 2] == 0.0
        @test stats.p_N[1, 1, 4] == 0.0
        @test abs(sum(stats.p_N[1, 1, :]) - 1.0) < 1e-12
    end

    @testset "fixed-N parity with the documented recipe" begin
        N_values = [0, 1, 2]
        ns_outputs = [DataFrame(iter=Int[], emax=Float64[]),
                      DataFrame(iter=collect(1:50), emax=zeros(50)),
                      DataFrame(iter=collect(1:60),
                                emax=collect(range(-0.010, -0.020, length=60)))]
        live_all = [Float64[], zeros(8), fill(-0.0201, 8)]
        T = 300.0
        mus = [mu_for(1.0, T), mu_for(3.0, T)]
        out = gc_thermodynamic_stats_fixed_N(ns_outputs, N_values, Vq, mass_ar,
            mus, [T]u"K"; n_walkers=8, live_emax=live_all)
        @test out.N_support == [0, 1, 2]
        @test size(out.p_N) == (2, 1, 3)
        β = 1.0 / (kb * T)
        for (k, mu) in enumerate(mus)
            w = out.log_Z_N[:, 1] .+ β * ustrip(u"eV", mu) .* N_values
            w .-= maximum(w)
            pw = exp.(w) ./ sum(exp.(w))
            for n_idx in 1:3
                @test isapprox(out.p_N[k, 1, n_idx], pw[n_idx]; rtol=1e-12, atol=1e-15)
            end
        end
    end
end
@testset "bounded-support normalization (truncated reference)" begin
    using Random
    # Calibration ledger:
    # - Bounded cross-route agreement (bounded igref reduction vs fixed-N
    #   stitching over N = 0..4, dilute LJ gas, K = 24, z0V = 6, n_max = 4,
    #   per-T mu grids solved from target zV in {4, 6}; igref seeds
    #   {1,2,3,62622} against fixed-N seeds {1001,1002,1003,63622}): max devs
    #   logXi 0.245, mean_N 0.530, N_eff floor observed 28.3 (binding on the
    #   T = 400 K reduction of the seed-3 calibration pair); gates ship at
    #   >= 3x (0.75 and 1.6) with an N_eff floor assertion at 10.
    # - Every other testset here is an exact identity (constructed ledgers or
    #   the same-input constant-shift contract) and needs no calibration.
    kb = 8.617333262e-5
    Vq = 1000.0u"Å^3"
    V = 1000.0
    mass_ar = 39.948u"u"
    lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(mass_ar, T * u"K"))
    mu_for(zV, T) = (kb * T * (log(zV / V) + 3 * log(lam(T)))) * u"eV"
    logfact(n) = n == 0 ? 0.0 : sum(log, 1:n)
    lpps = FreeBird.AnalysisTools._log_poisson_partial_sum

    # A compression ledger that ENCODES the truncated (conditional) Poisson
    # prior on 0:nmax exactly: shells carry the conditional masses of
    # N = 0..nmax-1 at energy zero, the tail walker carries N = nmax
    function trunc_poisson_ledger(z0V, nmax)
        logq = [n * log(z0V) - logfact(n) for n in 0:nmax]
        Q = sum(exp, logq)
        pi_n = exp.(logq) ./ Q
        X = reverse(cumsum(reverse(pi_n)))
        lc = [log(X[j+1] / X[j]) for j in 1:nmax]
        df = DataFrame(iter=collect(1:nmax), emax=zeros(nmax),
                       num_particles=collect(0:nmax-1), log_compression=lc)
        return df, [0.0], [nmax]
    end

    @testset "_log_poisson_partial_sum identities" begin
        @test lpps(3.0, 0) == 0.0
        @test abs(lpps(3.0, 60) - 3.0) < 1e-12
        @test abs(lpps(1.0, 1) - log(2.0)) < 1e-14
        @test lpps(3.0, 5) < lpps(3.0, 6) < lpps(3.0, 60)
        @test_throws ArgumentError lpps(3.0, -1)
    end

    @testset "exact closure of the bounded assembly" begin
        z0V = 3.0
        nmax = 5
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = trunc_poisson_ledger(z0V, nmax)
        T = 300.0
        zVs = [1.0, 3.0, 7.0]
        mus = [mu_for(zV, T) for zV in zVs]
        stats = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                                 live_emax=live_e, live_numbers=live_n,
                                                 n_max=nmax)
        for (i, zV) in enumerate(zVs)
            # the truncated model's exact grand potential: logXi = log sum_{N<=nmax} (zV)^N/N!
            @test abs(stats.logXi[i, 1] - lpps(zV, nmax)) < 1e-11 * max(lpps(zV, nmax), 1.0)
            # exact conditional-Poisson occupancy at the target activity
            logq = [n * log(zV) - logfact(n) for n in 0:nmax]
            q = exp.(logq) ./ sum(exp, logq)
            m1 = sum(q .* (0:nmax))
            @test abs(stats.mean_N[i, 1] - m1) < 1e-11 * max(m1, 1.0)
            for n in 0:nmax
                @test abs(stats.p_N[i, 1, n + 1] - q[n + 1]) < 1e-11
            end
            @test abs(sum(stats.p_N[i, 1, :]) - 1.0) < 1e-12
        end
    end

    @testset "constant-shift contract against the unbounded normalization" begin
        z0V = 3.0
        nmax = 5
        z0 = (z0V / V)u"Å^-3"
        df, live_e, live_n = trunc_poisson_ledger(z0V, nmax)
        T = 300.0
        mus = [mu_for(2.0, T), mu_for(5.0, T)]
        b = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                             live_emax=live_e, live_numbers=live_n,
                                             n_max=nmax)
        u = gc_thermodynamic_stats_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                             live_emax=live_e, live_numbers=live_n)
        # logXi differs by exactly the reference-mass swap; everything that the
        # reference mass cancels from is bitwise identical
        @test all(abs.((b.logXi .- u.logXi) .- (lpps(z0V, nmax) - z0V)) .< 1e-12)
        @test b.mean_N == u.mean_N
        @test b.var_N == u.var_N
        @test b.mean_U == u.mean_U
        @test b.p_N == u.p_N
        @test b.N_eff == u.N_eff
        # the ESS reduction needs no n_max: same ledger, same N_eff, bitwise
        ess = gc_effective_sample_size_ideal_ref(df, Vq, mass_ar, z0, mus, [T]u"K";
                                                 live_emax=live_e, live_numbers=live_n)
        @test ess == b.N_eff
    end

    @testset "bounded-run consistency validation" begin
        z0 = (3.0 / V)u"Å^-3"
        T = [300.0]u"K"
        mus = [mu_for(3.0, 300.0)]
        df, live_e, live_n = trunc_poisson_ledger(3.0, 5)
        # a ledger row above the claimed cap is rejected
        df_bad = copy(df)
        df_bad.num_particles[end] = 9
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df_bad, Vq, mass_ar, z0, mus, T;
            live_emax=live_e, live_numbers=live_n, n_max=5)
        # a live-tail count above the claimed cap is rejected
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, Vq, mass_ar, z0, mus, T;
            live_emax=live_e, live_numbers=[9], n_max=5)
        # a negative cap is rejected
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, Vq, mass_ar, z0, mus, T;
            live_emax=live_e, live_numbers=live_n, n_max=-1)
        # the typemax sentinel of an unbounded run is rejected with a clear
        # error instead of a runaway allocation
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, Vq, mass_ar, z0, mus, T;
            live_emax=live_e, live_numbers=live_n, n_max=typemax(Int64))
    end

    @testset "bounded cross-route agreement (seed 62622)" begin
        box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        pbc = (true, true, true)
        seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
        mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                               empty(position(seed_at, :)), empty(species(seed_at, :)),
                               empty(mass(seed_at, :)))
        V12 = 1728.0
        Vq12 = 1728.0u"Å^3"
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
        z0 = (6.0 / V12)u"Å^-3"
        n_cap = 4
        mu12(zV, T) = (kb * T * (log(zV / V12) + 3 * log(lam(T)))) * u"eV"
        Ts = [300.0, 400.0]
        mugrids = [[mu12(4.0, T), mu12(6.0, T)] for T in Ts]
        mksave(tag) = SaveEveryN(df_filename="_igbnd_$(tag).csv",
                                 wk_filename="_igbnd_$(tag).traj.extxyz",
                                 ls_filename="_igbnd_$(tag).ls.extxyz",
                                 n_traj=10^7, n_snap=10^7, n_info=10^7)
        clean(tag) = for f in ["_igbnd_$(tag).csv", "_igbnd_$(tag).traj.extxyz", "_igbnd_$(tag).ls.extxyz"]
            rm(f, force=true)
        end

        Random.seed!(62622)
        ls = GenericAtomWalkers([AtomWalker{1}(mkempty()) for _ in 1:24], lj)
        params = AtomisticIGRefGCNSParameters(mc_steps=60, reference_activity=z0,
                                              species=:Ar, allowed_fail_count=100,
                                              n_max=n_cap)
        df, lso, _ = ideal_gas_referenced_nested_sampling(
            ls, params, 300, MCAtomGrandCanonicalMoves(), mksave("ig"))
        clean("ig")
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]
        @test maximum(df.num_particles) <= n_cap
        @test maximum(live_n) <= n_cap

        Random.seed!(63622)
        dfs = DataFrame[]
        lives = Vector{Float64}[]
        for N in 0:n_cap
            if N == 0
                push!(dfs, DataFrame(iter=Int[], emax=Float64[]))
                push!(lives, zeros(8))
                continue
            end
            walkers = AtomWalker{1}[]
            for _ in 1:24
                coords = [[rand() * 12.0, rand() * 12.0, rand() * 12.0] for _ in 1:N]
                push!(walkers, AtomWalker{1}(FastSystem(atomic_system(
                    [:Ar => SVector{3}(c)u"Å" for c in coords], box, pbc))))
            end
            lsN = GenericAtomWalkers(walkers, lj)
            p = NestedSamplingParameters(mc_steps=60, initial_step_size=0.5, step_size=0.5,
                                         step_size_lo=0.01, step_size_up=2.0,
                                         allowed_fail_count=1000)
            dfN, lsoN, _ = nested_sampling(lsN, p, 150, MCRandomWalkClone(), mksave("fx$N"))
            clean("fx$N")
            push!(dfs, dfN)
            push!(lives, [ustrip(u"eV", w.energy) for w in lsoN.walkers])
        end

        for k in 1:2
            ig = gc_thermodynamic_stats_ideal_ref(df, Vq12, mass_ar, z0, mugrids[k], [Ts[k]]u"K";
                                                  live_emax=live_e, live_numbers=live_n,
                                                  n_max=n_cap)
            fx = gc_thermodynamic_stats_fixed_N(dfs, collect(0:n_cap), Vq12, mass_ar,
                                                mugrids[k], [Ts[k]]u"K";
                                                n_walkers=24, live_emax=lives)
            @test maximum(abs.(ig.logXi .- fx.logXi)) < 0.75
            @test maximum(abs.(ig.mean_N .- fx.mean_N)) < 1.6
            @test minimum(ig.N_eff) > 10.0
        end
    end
end
