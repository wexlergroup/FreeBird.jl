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
