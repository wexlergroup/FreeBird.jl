# Compression-convention coverage: the `compression` keyword of the iteration-based reductions
# (`ωᵢ`, the lattice ideal-gas-referenced stats and ESS, the fixed-N assembly, the microcanonical
# ladder) and the `compression` field of the serial atomistic parameter structs, with the `n_live`
# ledger column.
#
# Conventions. An ordinary cull from n live walkers leaves the fraction t ~ Beta(n, 1) of the enclosed
# prior volume. :geometric (the default) charges E[ln t] = -1/n, so ln X_i = -i/K and the shell above
# X_i is e^{-i/K}(e^{1/K} - 1); :mean charges ln E[t] = ln(n/(n+1)), the historical (K/(K+1))^i
# ladder, whose iteration-based shell is X_i - X_{i+1} and needs omega0 = (K+1)/K to close with the
# live tail. Plateau-tie evictions charge ln((n-1)/n) under both.
#
# Calibration ledger (gates at >= 3x the largest three-seed deviation, per statistic; the seed shipped
# in the testset is the first calibration seed, so its statistics reproduce digit for digit; run with
# FB_CALIBRATE=1 to print all three):
# - Synthetic ledger at depth (K = 50, n = 7000, D = 100, s = 4, R = 200 replicates; seeds
#   424601/424602/424603): mean convention bias +0.8653, +0.8652, +0.8964 nats at z = 11.06, 9.86, 9.54
#   (gate: z > 3); geometric bias +0.0466, +0.0454, +0.0775 at z = 0.59, 0.51, 0.82 (gate: |z| < 2.5);
#   shift between the conventions 0.8187, 0.8198, 0.8189 against the third-order prediction 0.8189
#   (the shells' prefactor ratio plus (j_dom - 1)(1/K - ln(1 + 1/K)); deviations 0.0003, 0.0008, 0.0000;
#   gate: |deviation| < 0.003); positive fraction of the geometric deviations 0.535,
#   0.505, 0.540 (gate: (0.38, 0.62)); sd ratio geometric/mean 1.0095, 1.0089, 1.0088 (the spreads
#   agree to first order; loose band (0.8, 1.25), a report).
# - The lattice forwarding testset gates exact identities only (rtol 1e-10), no statistics.
const CC_GATE_ZM = 3.0             # mean-convention excess in standard errors (calibrated z 9.54 to 11.06)
const CC_GATE_ZG = 2.5             # |geometric bias| band in standard errors (calibrated |z| 0.51 to 0.82)
const CC_GATE_SHIFT = 0.003        # nats, |mean shift - the third-order prediction| (calibrated 0.0000 to 0.0008)
const CC_GATE_FRAC = (0.38, 0.62)  # positive fraction of the geometric deviations (calibrated 0.505 to 0.540)
const CC_GATE_FWD = 1e-12          # rtol of the exact per-shell ratio identity on the absolute quantities (no statistics)
const CC_CALIBRATE = get(ENV, "FB_CALIBRATE", "") == "1"

@testset "compression convention" begin
    using Random
    using Statistics
    using DataFrames
    using Unitful

    # ---------------------------------------------------------------- 1. shell identities
    @testset "shell identities under both conventions" begin
        for K in (4, 32, 250), n in (1, 10, 1000)
            sg = ωᵢ(collect(1:n), K)                              # default: geometric
            sm = ωᵢ(collect(1:n), K; compression=:mean)
            @test sg == ωᵢ(collect(1:n), K; compression=:geometric)
            # geometric shells are X_{i-1} - X_i with ln X_i = -i/K: they sum to 1 - X_n at omega0 = 1
            @test isapprox(sum(sg), 1 - exp(-n / K); rtol=1e-12)
            # mean shells are X_i - X_{i+1}: they sum to r(1 - r^n), so omega0 = (K+1)/K closes them
            r = K / (K + 1)
            @test isapprox(sum(sm), r * (1 - r^n); rtol=1e-12)
            @test isapprox(sum(ωᵢ(collect(1:n), K; ω0=(K + 1) / K, compression=:mean)) + r^n, 1.0; rtol=1e-12)
            # the per-shell ratio and the cumulative-mass ratio between the conventions
            for i in (1, min(10, n), n)
                @test isapprox(sg[i] / sm[i], exp(-i / K) * expm1(1 / K) / ((1 / (K + 1)) * r^i); rtol=1e-12)
            end
            # cumulative masses: ln X_n(mean) - ln X_n(geometric) = n (ln r + 1/K) > 0, read off the
            # closed shells where the tail is resolvable (X_n > 1e-6)
            if r^n > 1e-6
                smc = ωᵢ(collect(1:n), K; ω0=(K + 1) / K, compression=:mean)
                @test isapprox(log(1 - sum(smc)) - log(1 - sum(sg)), n * (log(r) + 1 / K); rtol=1e-6)
            end
            @test all(diff(sg) .< 0) && all(sg .> 0)
        end
        # n_cull > 1: the geometric branch charges -C/K per iteration
        @test isapprox(sum(ωᵢ(collect(1:50), 20; n_cull=3)), 1 - exp(-50 * 3 / 20); rtol=1e-12)
        # the mean branch is the historical expression, bit for bit
        @test ωᵢ([1, 2, 3], 4; compression=:mean) == 1.0 * (1 / (4 + 1)) * (4 / (4 + 1)) .^ [1, 2, 3]
        # edge cases behave alike
        @test ωᵢ(Int[], 4) == Float64[] && ωᵢ(Int[], 4; compression=:mean) == Float64[]
        @test_throws MethodError ωᵢ([1.5], 4; compression=:geometric)
        @test_throws ArgumentError ωᵢ([1], 4; compression=:harmonic)
        @test_throws ArgumentError NestedSamplingParameters(compression=:harmonic)
        @test_throws ArgumentError AtomisticIGRefGCNSParameters(compression=:harmonic)
        @test NestedSamplingParameters().compression === :geometric
        @test AtomisticIGRefGCNSParameters().compression === :geometric
        @test NestedSamplingParameters(compression=:mean).compression === :mean
    end

    # ---------------------------------------------------------------- 2. the bias over runs, exactly solvable
    @testset "synthetic ledger at depth: the mean convention's ln Z excess (seeded, calibrated)" begin
        # Shrinkage factors t_j ~ Beta(K, 1) drawn as exp(ln u / K) from raw uniforms (no Distributions
        # sampler), X_j = prod t; a likelihood peaked at depth D in ln X, L(X) = exp(-(ln X + D)^2 / (2 s^2)),
        # whose evidence is closed: Z = e^{s^2/2 - D} s sqrt(2 pi) Phi((D - s^2)/s). Each replicate reduces
        # the SAME dead points under both conventions; the estimator is unbiased for the geometric masses
        # to first order and carries about +(D - s^2)/(2K) under the mean masses (the weight of e^u L(u)
        # peaks at u = -D + s^2, so the dominant index is j ~ K (D - s^2)).
        Φ(x) = 0.5 * (1 + erf_approx(x / sqrt(2)))
        # erf is not in Base and SpecialFunctions is not a test dependency: the Numerical Recipes erfc
        # Chebyshev fit (accurate to about 1e-7) is more than the test needs at (D - s^2)/s = 21, where
        # log(Phi) evaluates to exactly 0.0.
        function erf_approx(x)
            z = abs(x); t = 1 / (1 + 0.5 * z)
            r = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
                t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
                t * (-0.82215223 + t * 0.17087277)))))))))
            return x >= 0 ? 1 - r : r - 1
        end
        function synthetic(seed; K=50, n=7000, D=100.0, s=4.0, R=200)
            rng = MersenneTwister(seed)
            lnZ_true = s^2 / 2 - D + log(s * sqrt(2π)) + log(Φ((D - s^2) / s))
            devm = Float64[]; devg = Float64[]
            for _ in 1:R
                lnX = cumsum(log.(rand(rng, n)) ./ K)               # the true ln X_j
                lnL = -(lnX .+ D) .^ 2 ./ (2 * s^2)                   # ln L at the dead points
                for (conv, ω0) in ((:mean, (K + 1) / K), (:geometric, 1.0))
                    w = ωᵢ(collect(1:n), K; ω0=ω0, compression=conv)
                    Xn = conv === :mean ? (K / (K + 1))^n : exp(-n / K)
                    lt = vcat(log.(w) .+ lnL, fill(log(Xn) - (lnX[end] + D)^2 / (2 * s^2), 1))  # tail at X_n
                    m = maximum(lt)
                    push!(conv === :mean ? devm : devg, m + log(sum(exp.(lt .- m))) - lnZ_true)
                end
            end
            # the shift between the conventions to third order: the shells' prefactor ratio
            # -log((K+1)(1 - e^{-1/K})) plus (j_dom - 1) times the per-cull difference 1/K - log(1 + 1/K)
            # at the dominant index j_dom = K (D - s^2) (the crude (D - s^2)/(2K) is high by about 1/(2K))
            predicted = -log((K + 1) * (-expm1(-1 / K))) + (K * (D - s^2) - 1) * (1 / K - log1p(1 / K))
            return (mean_bias=mean(devm), se_m=std(devm) / sqrt(R), geo_bias=mean(devg), se_g=std(devg) / sqrt(R),
                    frac_pos_g=count(>(0), devg) / R, frac_pos_m=count(>(0), devm) / R,
                    shift=mean(devm .- devg), sd_ratio=std(devg) / std(devm), predicted=predicted)
        end
        if CC_CALIBRATE
            for seed in (424601, 424602, 424603)
                c = synthetic(seed)
                println("CALIBRATION synthetic seed $seed: mean_bias $(c.mean_bias) se_m $(c.se_m) z_m $(c.mean_bias / c.se_m) geo_bias $(c.geo_bias) se_g $(c.se_g) z_g $(c.geo_bias / c.se_g) frac_pos_g $(c.frac_pos_g) frac_pos_m $(c.frac_pos_m) shift $(c.shift) predicted $(c.predicted) sd_ratio $(c.sd_ratio)")
            end
        end
        c = synthetic(424601)
        # the mean convention sits above the truth by the third-order prediction (about (D - s^2)/(2K))
        # within CC_GATE_SHIFT
        @test c.mean_bias > 0
        @test c.mean_bias > CC_GATE_ZM * c.se_m
        @test abs(c.shift - c.predicted) < CC_GATE_SHIFT
        # the geometric convention is sign-balanced: its mean within CC_GATE_ZG standard errors of zero
        # and its positive fraction inside the calibrated band
        @test abs(c.geo_bias) < CC_GATE_ZG * c.se_g
        @test CC_GATE_FRAC[1] < c.frac_pos_g < CC_GATE_FRAC[2]
        @test c.frac_pos_m > c.frac_pos_g
        # the run-to-run spreads agree to first order (the two reductions differ by a near-constant
        # tilt); reported, not gated beyond a loose band
        @test 0.8 < c.sd_ratio < 1.25
    end

    # ---------------------------------------------------------------- 3. the sampler charge and the n_live column
    @testset "serial atomistic steps: charge and n_live under both conventions (tie fixture)" begin
        cc_L = 12.0
        cc_box = [[cc_L * u"Å", 0u"Å", 0u"Å"], [0u"Å", cc_L * u"Å", 0u"Å"], [0u"Å", 0u"Å", cc_L * u"Å"]]
        cc_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=3.0, shift=true)
        cc_pair(r) = AtomWalker(FastSystem(periodic_system(
            [:Ar => [0.3, 0.5, 0.5], :Ar => [0.3 + r / cc_L, 0.5, 0.5]], cc_box, fractional=true)))
        cc_save(tag) = SaveEveryN(df_filename="_cc_$(tag).csv", wk_filename="_cc_$(tag).w.extxyz",
                                  ls_filename="_cc_$(tag).l.extxyz", n_traj=10^8, n_snap=10^8, n_info=10^8)
        cc_rm(tag) = for f in ("_cc_$(tag).csv", "_cc_$(tag).w.extxyz", "_cc_$(tag).l.extxyz")
            rm(f, force=true)
        end
        # three bit-exact duplicate dimers above two deeper ones: a 3-tie block, then ordinary culls
        function tie_run(compression, tag)
            Random.seed!(76543)
            dup = cc_pair(2.90)
            ls = LJAtomWalkers([deepcopy(dup), deepcopy(dup), deepcopy(dup), cc_pair(2.82), cc_pair(2.84)], cc_lj)
            p = NestedSamplingParameters(mc_steps=5, initial_step_size=0.5, step_size=0.5,
                step_size_lo=0.01, step_size_up=2.0, accept_range=(0.25, 0.75),
                allowed_fail_count=1000, energy_perturbation=1e-12, compression=compression)
            df, fls, _ = nested_sampling(ls, p, 8, MCGalileanWalk(n_refresh=4), cc_save(tag))
            cc_rm(tag)
            return df, fls
        end
        dfg, flg = tie_run(:geometric, "g")
        dfm, flm = tie_run(:mean, "m")
        @test names(dfg) == ["iter", "emax", "log_compression", "n_live"]
        @test names(dfm) == names(dfg)
        # identical trajectories: the charge never enters a walk or an acceptance decision
        @test dfg.iter == dfm.iter && dfg.emax == dfm.emax && dfg.n_live == dfm.n_live
        @test [ustrip(u"eV", w.energy) for w in flg.walkers] == [ustrip(u"eV", w.energy) for w in flm.walkers]
        # the tie block: three evictions from 5, 4, 3 live walkers, charged (n-1)/n under both
        @test dfg.n_live[1:3] == [5, 4, 3]
        @test dfg.log_compression[1:3] ≈ [log(4 / 5), log(3 / 4), log(2 / 3)] atol = 1e-14
        @test dfm.log_compression[1:3] == dfg.log_compression[1:3]
        # ordinary culls after the refill: -1/n_live versus log(n_live/(n_live + 1)), n_live recorded
        ord = 4:nrow(dfg)
        @test !isempty(ord)
        @test all(dfg.log_compression[ord] .== -1.0 ./ dfg.n_live[ord])
        @test all(dfm.log_compression[ord] .== log.(dfm.n_live[ord] ./ (dfm.n_live[ord] .+ 1)))
        @test all(dfg.n_live[ord] .<= 5) && all(dfg.n_live[ord] .>= 2)
        # both ledgers reduce through the same ledger method and both telescope with their tails
        for df in (dfg, dfm)
            w = ωᵢ(Vector{Float64}(df.log_compression))
            @test isapprox(sum(w) + exp(sum(df.log_compression)), 1.0; rtol=1e-12)
        end
        # exact conversion through n_live: the mean ledger maps onto the geometric one row by row
        conv = [dfm.log_compression[i] == log(dfm.n_live[i] / (dfm.n_live[i] + 1)) ? -1.0 / dfm.n_live[i] :
                dfm.log_compression[i] for i in 1:nrow(dfm)]
        @test conv == dfg.log_compression
        # the record_move_rates schema declares the column eagerly, right after log_compression
        Random.seed!(76545)
        ls3 = LJAtomWalkers([cc_pair(2.80 + 0.02k) for k in 1:4], cc_lj)
        p3 = NestedSamplingParameters(mc_steps=5, initial_step_size=0.5, step_size=0.5,
            step_size_lo=0.01, step_size_up=2.0, accept_range=(0.25, 0.75),
            allowed_fail_count=1000, energy_perturbation=1e-12)
        df3, _, _ = nested_sampling(ls3, p3, 3, MCRandomWalkClone(), cc_save("r"); record_move_rates=true)
        cc_rm("r")
        @test names(df3)[1:4] == ["iter", "emax", "log_compression", "n_live"]
        @test all(df3.log_compression .== -1.0 ./ df3.n_live)
        # the column is reserved against observable-name collisions
        @test :n_live in FreeBird.SamplingSchemes._RESERVED_LEDGER_COLUMNS
    end

    # ---------------------------------------------------------------- 4. forwarding on the lattice fixture
    @testset "lattice ideal-gas-referenced stats forward the keyword (exact identities)" begin
        # the 4 x 4 IDEAL lattice gas (zero couplings) at z0 = 1: every configuration has E = 0, so at the
        # reference point the reduction is the mass law alone and ln Xi = M ln 2 under both idioms;
        # away from it the two conventions differ by the exact weighted per-shell ratio, which the
        # reduction must reproduce from omega_i itself (the forwarding statement, no statistics)
        M = 16
        lat = MLattice{1,SquareLattice}(lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 1), periodicity=(true, true, false), cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:M]], adsorptions=:full)
        ham = GenericLatticeHamiltonian(0.0, [0.0, 0.0], u"eV")
        kb = 8.617333262e-5
        Random.seed!(424611)
        K = 16
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:K]
        ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=50, reference_fugacity=1.0)
        save = SaveEveryN("_cc_lat.csv", "_cc_lat.traj", "_cc_lat.ls", 10^6, 10^6, 10^6)
        df, fls, _ = ideal_gas_referenced_nested_sampling(ls, params, 300, MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3), save)
        foreach(f -> rm(f, force=true), ("_cc_lat.csv", "_cc_lat.traj", "_cc_lat.ls"))
        live_E = [w.energy.val for w in fls.walkers]
        live_N = [sum(w.configuration.components[1]) for w in fls.walkers]
        n = nrow(df)
        @test n > 50
        μ = 0.0; T = 300.0
        sg = gc_thermodynamic_stats_ideal_ref(df, M, 1.0, [μ], [T], K; live_emax=live_E, live_numbers=live_N)
        sd = gc_thermodynamic_stats_ideal_ref(df, M, 1.0, [μ], [T], K; live_emax=live_E, live_numbers=live_N, compression=:geometric)
        sm = gc_thermodynamic_stats_ideal_ref(df, M, 1.0, [μ], [T], K; ω0=(K + 1) / K, live_emax=live_E,
                                              live_numbers=live_N, compression=:mean)
        @test sg.logXi == sd.logXi && sg.mean_N == sd.mean_N && sg.N_eff == sd.N_eff
        # the mass law at the reference point under both idioms
        @test isapprox(sg.logXi[1, 1], M * log(2.0); atol=1e-10)
        @test isapprox(sm.logXi[1, 1], M * log(2.0); atol=1e-10)
        # away from the reference point: ln Xi(mean) - ln Xi(geometric) equals the log of the weighted
        # per-shell ratio, with the weights the geometric reduction assigns and the ratios from omega_i
        μ2 = 0.03
        sg2 = gc_thermodynamic_stats_ideal_ref(df, M, 1.0, [μ2], [T], K; live_emax=live_E, live_numbers=live_N)
        sm2 = gc_thermodynamic_stats_ideal_ref(df, M, 1.0, [μ2], [T], K; ω0=(K + 1) / K, live_emax=live_E,
                                               live_numbers=live_N, compression=:mean)
        β = 1 / (kb * T)
        wg = ωᵢ(collect(1:n), K)
        wm = ωᵢ(collect(1:n), K; ω0=(K + 1) / K, compression=:mean)
        lw = vcat(log.(wg), fill(-n / K - log(K), K))
        ratio = vcat(wm ./ wg, fill((K / (K + 1))^n / exp(-n / K), K))
        Ns = vcat(Float64.(df.num_particles), Float64.(live_N))
        lt = lw .+ β * μ2 .* Ns
        w = exp.(lt .- maximum(lt)); w ./= sum(w)
        # (the sign of the difference depends on where the weight sits: the tilt j/(2K^2) competes
        # with the shells' constant prefactor ratio 1 - 1/(2K) on a shallow ladder; the identity is
        # gated on the absolute quantities, since the difference itself is a small number formed from
        # two O(10) log-sum-exp reductions along different code paths)
        @test isapprox(sm2.logXi[1, 1], sg2.logXi[1, 1] + log(sum(w .* ratio)); rtol=CC_GATE_FWD)
        # the ESS reduction forwards the keyword: the geometric default equals the explicit keyword and
        # differs from the mean idiom away from the reference point
        eg = gc_effective_sample_size_ideal_ref(df, M, 1.0, [μ2], [T], K; live_emax=live_E, live_numbers=live_N)
        ed = gc_effective_sample_size_ideal_ref(df, M, 1.0, [μ2], [T], K; live_emax=live_E, live_numbers=live_N, compression=:geometric)
        em = gc_effective_sample_size_ideal_ref(df, M, 1.0, [μ2], [T], K; ω0=(K + 1) / K, live_emax=live_E,
                                                live_numbers=live_N, compression=:mean)
        @test eg == ed
        @test eg != em
        # the fixed-N assembly: the geometric default with omega0 = 1 gives a sector mass of exactly one
        # (a single sector with a flat ladder at E = 0 carries ln Z_N = ln C(M, N)), the mean idiom the
        # historical value
        dfs = [DataFrame(iter=Int[], emax=Float64[]), DataFrame(iter=collect(1:30), emax=zeros(30))]
        fx = gc_thermodynamic_stats_fixed_N(dfs, [0, 1], M, [0.0u"eV"], [300.0u"K"]; n_walkers=10,
                                            live_emax=[[0.0], zeros(10)])
        fxm = gc_thermodynamic_stats_fixed_N(dfs, [0, 1], M, [0.0u"eV"], [300.0u"K"]; n_walkers=10,
                                             ω0=11 / 10, live_emax=[[0.0], zeros(10)], compression=:mean)
        @test isapprox(fx.log_Z_N[2, 1], log(M); rtol=1e-12)
        @test isapprox(fxm.log_Z_N[2, 1], log(M) + log(1 + (10 / 11)^30 / 10); rtol=1e-12)
        # cv(df, ...) and the omega-column gc_thermodynamic_stats(df, ...) forward the keyword: the
        # mean-convention weights of a flat ledger are the hand-built (1/(K+1)) r^i, and the geometric
        # default differs from them
        dfc = DataFrame(iter=collect(1:40), emax=collect(range(1.0, 0.0; length=40)))
        cvg = cv(dfc, [30.0, 40.0], 0, 10)
        cvm = cv(dfc, [30.0, 40.0], 0, 10; compression=:mean)
        @test cvg == cv(dfc, [30.0, 40.0], 0, 10; compression=:geometric)
        @test cvm != cvg
        wm_hand = (1 / 11) .* (10 / 11) .^ dfc.iter
        @test cvm ≈ [cv(b, wm_hand, dfc.emax .- minimum(dfc.emax), 0) for b in [30.0, 40.0]] rtol = 1e-12
        dfo = DataFrame(iter=collect(1:40), omega=collect(range(1.0, 0.0; length=40)),
                        energy=collect(range(1.0, 0.0; length=40)), num_particles=collect(1:40))
        gg = gc_thermodynamic_stats(dfo, [30.0, 40.0], 10, 0.0)
        gm = gc_thermodynamic_stats(dfo, [30.0, 40.0], 10, 0.0; compression=:mean)
        @test gg == gc_thermodynamic_stats(dfo, [30.0, 40.0], 10, 0.0; compression=:geometric)
        @test gm[3] != gg[3]
        # the microcanonical ladder slope follows the keyword: the volume entropies are proportional
        dfl = DataFrame(iter=collect(1:400), emax=collect(range(10.0, 0.0; length=400)))
        Eg, Sg = microcanonical_entropy(dfl, 50; kind=:volume)
        Emn, Smn = microcanonical_entropy(dfl, 50; kind=:volume, compression=:mean)
        @test Eg == Emn
        @test isapprox(maximum(abs.((Sg .- Sg[end]) .- (Smn .- Smn[end]) .* ((-1 / 50) / log(50 / 51)))), 0.0; atol=1e-9)
    end
end
