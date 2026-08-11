# Regression coverage for plateau-aware tie eviction in the serial atomistic
# nested-sampling steps (issue #181, option 1) and the log-compression method of
# ωᵢ. Exact energy ties at the ceiling arise generically from truncated pair
# potentials over configuration spaces with vacuum (lj_energy is exactly
# 0.0 eV beyond the cutoff); the fix evicts tied walkers one by one without
# replacement, compressing by (n_live - 1)/n_live with the shrinking live count
# (Fowlie, Handley and Su, Mon. Not. R. Astron. Soc. 503, 1199 (2021)), records
# the per-cull log-compression in a lazily added ledger column, and refills the
# live set only once the plateau is exhausted.
#
# Fixture: the reduced-scale port of issue #181's inline reproducer, sites
# verbatim: a frozen 13-atom cuboctahedral LJ cluster (eps = 0.01 eV,
# sigma = 2.5 A, cutoff 3 sigma shifted) centered in an all-false-periodicity
# cubic box padded by two cutoff radii of vacuum, with N ideal adsorbates
# (adsorbate-adsorbate epsilon = 0 through CompositeParameterSets),
# LJSurfaceWalkers, K = 48, MCRandomWalkClone. The exact per-N reference is
# N*ln(z1) with z1 a midpoint-quadrature single-particle integral over the same
# box, restricted to the walker-initialization prior E < 1e9 eV (the "m" fold
# below divides by the admissible-point count, matching uniform_walker's
# rejection). The surface's energy_frozen_part is set to 0 eV by hand, so walker
# energies carry the adsorbate part only and no beta*E_frozen fold enters the
# closure.
#
# ================= calibration ledger (record, do not retype) =================
# All statistical gates below ship at >= 3x the maximum multi-seed deviation of
# their own statistic (per stat, per N sector), calibrated on the shipped seeds
# themselves at T = 100 K. Seeds are the first-scanned blocks 11-13 (N = 1),
# 21-23 (N = 2), 41-43 (control); they were not selected on their deviations.
# Digit-for-digit deviations (reproduced by running this file with
# FREEBIRD_NS_PLATEAU_TIES_CALIBRATE=1; gates evaluate after all RNG use, so
# comment/tolerance edits never invalidate the calibration):
#
#   N=1 (mc_steps=400, n_steps=400, K=48), dev_ev = lnZ - ln z1 at 100 K,
#   dev_mass = ln(recorded plateau mass / f):
#     seed 11: n_tie=40  dev_ev_new=-0.0013862869792976226  dev_ev_leg=+0.08734514794562426
#              dev_mass_new=-0.08575642023567279  dev_mass_leg=-0.48028481527974104
#     seed 12: n_tie=44  dev_ev_new=-0.02480029212362888  dev_ev_leg=+0.09662588733381072
#              dev_mass_new=+0.009553759568652065  dev_mass_leg=-0.42033429427638
#     seed 13: n_tie=43  dev_ev_new=-0.07190702381754824  dev_ev_leg=-0.04934557228207022
#              dev_mass_new=-0.03405504585878227  dev_mass_leg=-0.45515436372078055
#              (seed 13 culls one positive-energy walker before the block, so
#              its block statistic carries one ordinary-cull X factor)
#     GATE_EV_N1   = 0.2158 >= 3 x 0.07190702 = 0.21572107 (max |dev_ev_new|, seed 13)
#     GATE_MASS_N1 = 0.2573 >= 3 x 0.08575642 = 0.25726926 (max |dev_mass_new|, seed 11)
#     contrast margins on the mechanism ledger (seed 11):
#       |dev_mass_leg| = 0.48028 = 1.87 x GATE_MASS_N1 (fails low: plateau mass
#       under-assigned); sub-plateau mass overshoot ln(S_leg/S_true) = +1.56
#       (fails high; structurally >= +1.37 for any n_tie <= K at this f,
#       asserted > 1.0).
#
#   N=2 (mc_steps=400, n_steps=600, K=48):
#     seed 21: n_tie=38  dev_ev_new=-0.08329687661812356  dev_ev_leg=+0.020386112867883957
#              dev_mass_new=-0.16420030128135582
#     seed 22: n_tie=38  dev_ev_new=-0.04519088682375805  dev_ev_leg=+0.03614436441601944
#              dev_mass_new=-0.040484578064941684  (run ends inside a later
#              micro-tie block: live=45, plateau_refill_target=48; the end-state
#              contract assert covers both endings)
#     seed 23: n_tie=42  dev_ev_new=-0.08820427561484093  dev_ev_leg=+0.02154990570166611
#              dev_mass_new=+0.03897959328930537
#     GATE_EV_N2 = 0.2647 >= 3 x 0.08820428 = 0.26461283 (max |dev_ev_new|, seed 23)
#     The N=2 block-mass statistic is recorded (printed under the calibration
#     flag) but not gated: several positive-energy pre-block culls can shift
#     the block's starting prior mass (seed 21 drew -0.164), so an honest
#     >= 3x gate would be too loose to mean anything at N=2, and the
#     accounting it would probe is already locked at N=1.
#
#   control (cutoff=Inf, N=1, mc_steps=400, n_steps=400):
#     seed 41: dev_ev_new=-0.0068748326320041825
#     seed 42: dev_ev_new=+0.12225112083863904
#     seed 43: dev_ev_new=-0.05334596167653882
#     GATE_EV_CTRL = 0.3668 >= 3 x 0.12225112 = 0.36675336 (max |dev|, seed 42)
#
# Choice of closure sectors: N in {1, 2} rather than the issue's {1, 4}; N = 4
# needs a ladder roughly 4x deeper for the same tail control and exercises no
# accounting path beyond N = 2 (which already produces secondary micro-tie
# blocks from clone energy-residue degeneracy). Gate temperature 100 K: the
# 400-row ladders close tightly there (Boltzmann weight below the recorded
# ladder is negligible and the live-set tail term covers the remainder). At
# T <= 60 K the same ledgers under-estimate the evidence by up to ~1.5 nats
# REGARDLESS of weighting: after the plateau is exhausted the refill clones the
# few sub-plateau survivors (~K(1-f) ~ 4 walkers here), and the subsequent E(X)
# mapping is ancestor-correlated for ~1-2 e-folds of compression. That is a
# sampler-variance property of the fixture scale, not an accounting defect, and
# it is why the legacy-vs-corrected contrast is asserted on the recorded
# plateau prior mass (whose only noise is the binomial live-set draw, +-0.09
# here) rather than on the low-T evidence, where the two effects are
# comparable in size (verified over 20-seed scans at T = 40-300 K).
# ==============================================================================

using Random

@testset "Atomistic NS plateau ties (issue #181)" begin

    calib = get(ENV, "FREEBIRD_NS_PLATEAU_TIES_CALIBRATE", "0") == "1"

    # ================= item 2: deterministic weight identities =================
    # Synthetic log-compression vectors only; no sampling, no RNG.
    @testset "log-compression ωᵢ: deterministic weight identities" begin
        K = 32
        n_pre, n_tie, n_post = 7, 12, 9
        # n_pre ordinary culls, a hand-placed n_tie-row eviction block with the
        # shrinking live count, n_post ordinary culls.
        logc = vcat(fill(log(K / (K + 1)), n_pre),
                    [log((K - 1 - j) / (K - j)) for j in 0:(n_tie - 1)],
                    fill(log(K / (K + 1)), n_post))
        shells = ωᵢ(logc)
        @test all(shells .> 0)

        # (a) the tie block reproduces the closed-form block factor
        # (K - n_tie)/K: telescoping product of (n-1)/n for n = K..K-n_tie+1.
        # Both routes measured exactly 0.0 at calibration; gates 1e-14.
        X_pre = 1 - sum(shells[1:n_pre])
        X_post = 1 - sum(shells[1:n_pre+n_tie])
        @test abs(X_post / X_pre - (K - n_tie) / K) <= 1e-14
        @test abs(exp(sum(logc[n_pre+1:n_pre+n_tie])) - (K - n_tie) / K) <= 1e-14
        # hand-built 7-tie block from 16 live walkers: factor 9/16 (measured 0.0)
        K7, n7 = 16, 7
        logc7 = [log((K7 - 1 - j) / (K7 - j)) for j in 0:(n7 - 1)]
        @test abs(exp(sum(logc7)) - (K7 - n7) / K7) <= 1e-14

        # (b) tie-free uniform column: the new method with its ω0 = 1.0 default
        # matches the legacy iteration-based weights called as
        # ωᵢ(1:n, K; ω0=(K+1)/K). Calibrated max per-element rel dev 1.705e-14.
        n = 250
        shells_u = ωᵢ(fill(log(K / (K + 1)), n))
        shells_l = ωᵢ(collect(1:n), K; ω0=(K + 1) / K)
        @test all(shells_u .> 0)
        @test all(abs.(shells_u .- shells_l) .<= 1e-12 .* shells_l)

        # (c) legacy-method pins on fixed inputs: the iteration-based code path
        # is untouched by the fix. Closed form ω0*(1/(K+1))*(K/(K+1))^i at
        # K = 48, ω0 = 49/48, i in (1, 5, 20): exactly (1/49)*(48/49)^(i-1);
        # digits below are the shipped code path's evaluation (repr, dev SHA
        # cfa19b8 semantics).
        w_leg = ωᵢ([1, 5, 20], 48; ω0=49 / 48)
        @test isapprox(w_leg[1], 0.020408163265306117; rtol=1e-13)
        @test isapprox(w_leg[2], 0.018792499586397383; rtol=1e-13)
        @test isapprox(w_leg[3], 0.013793100784800585; rtol=1e-13)

        # (d) shells all positive (above) and total mass + tail identity:
        # sum(shells) + exp(sum(log_compression)) = 1 (measured 0.0).
        @test abs(sum(shells) + exp(sum(logc)) - 1.0) <= 1e-12
    end

    # ================= shared fixture: issue #181 reproducer =================
    kB = 8.617333262e-5              # eV K^-1

    sig  = 2.5                       # Angstrom
    eps0 = 0.01                      # eV
    rcut = 3.0                       # cutoff, units of sigma
    d_nn = 2.0^(1 / 6) * sig         # LJ pair-minimum distance
    L    = 2 * d_nn + 4 * rcut * sig # cluster diameter + two cutoff radii per side
    half = L / 2

    # 13-atom cuboctahedral cluster sites, verbatim from the issue #181 body.
    fcc_shell = [(1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),
                 (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),
                 (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1)]
    sites = vcat([(half, half, half)],
                 [(half + p[1] * d_nn / sqrt(2), half + p[2] * d_nn / sqrt(2),
                   half + p[3] * d_nn / sqrt(2)) for p in fcc_shell])

    lj_field(cut) = LJParameters(epsilon=eps0, sigma=sig, cutoff=cut, shift=true)
    lj_ideal(cut) = LJParameters(epsilon=0.0, sigma=sig, cutoff=cut, shift=false)
    # [adsorbate-adsorbate (ideal), adsorbate-surface, surface-surface]
    plateau_cps(cut) = CompositeParameterSets(2, [lj_ideal(cut), lj_field(cut), lj_field(cut)])

    ties_box = [[L, 0.0, 0.0]u"Å", [0.0, L, 0.0]u"Å", [0.0, 0.0, L]u"Å"]
    ties_pbc = (false, false, false)

    surface_system() = FastSystem(atomic_system(
        [:Ar => [x, y, z]u"Å" for (x, y, z) in sites], ties_box, ties_pbc))

    # single-adsorbate field energy at a point, through the library's own pair fn
    function field_energy(x, y, z, lj)
        e = 0.0u"eV"
        for (sx, sy, sz) in sites
            r = sqrt((x - sx)^2 + (y - sy)^2 + (z - sz)^2)u"Å"
            e += lj_energy(r, lj)
        end
        return ustrip(u"eV", e)
    end

    E_THRESH = 1.0e9                 # walker-initialization admissibility, eV
    NGRID = 144                      # midpoint quadrature, as in the issue body

    function field_grid(lj)
        h = L / NGRID
        E = Vector{Float64}(undef, NGRID^3)
        k = 0
        for iz in 1:NGRID, iy in 1:NGRID, ix in 1:NGRID
            k += 1
            E[k] = field_energy((ix - 0.5) * h, (iy - 0.5) * h, (iz - 0.5) * h, lj)
        end
        return E
    end

    # the m fold: mean over the admissible prior E < E_THRESH, matching the
    # uniform_walker rejection below (issue reproducer identity, ported)
    function ref_ln_z1(E, beta)
        z = 0.0
        m = 0
        for e in E
            e < E_THRESH || continue
            m += 1
            z += exp(-beta * e)
        end
        return log(z / m)
    end

    plateau_fraction(E) = count(iszero, (e for e in E if e < E_THRESH)) /
                          count(<(E_THRESH), E)

    function uniform_walker(N, ljs, surface)
        while true
            coor = [:Ar => [rand() * L, rand() * L, rand() * L]u"Å" for _ in 1:N]
            w = AtomWalker(FastSystem(atomic_system(coor, ties_box, ties_pbc)))
            E = ustrip(u"eV", interacting_energy(w.configuration, ljs, w.list_num_par,
                                                 w.frozen, surface.configuration))
            isfinite(E) && E < E_THRESH && return w
        end
    end

    function run_plateau_ns(N, ljs, surface, seed; mc_steps=400, n_steps=400, K=48)
        Random.seed!(seed)   # the parameters' random_seed field is not consumed
        walkers = [uniform_walker(N, ljs, surface) for _ in 1:K]
        ls = LJSurfaceWalkers(walkers, ljs, surface)
        params = NestedSamplingParameters(mc_steps=mc_steps, step_size=1.0,
                                          step_size_up=3.0, allowed_fail_count=100_000)
        save = SaveEveryN(df_filename="_test_ties_df.csv",
                          wk_filename="_test_ties.traj.extxyz",
                          ls_filename="_test_ties.ls.extxyz",
                          n_traj=10_000_000, n_snap=10_000_000, n_info=10_000_000)
        df, ls_out, params_out = nested_sampling(ls, params, n_steps,
                                                 MCRandomWalkClone(), save)
        rm("_test_ties_df.csv", force=true)
        rm("_test_ties.traj.extxyz", force=true)
        rm("_test_ties.ls.extxyz", force=true)
        live = [ustrip(u"eV", w.energy) for w in ls_out.walkers]
        return df, live, params_out
    end

    # per-N evidence, new accounting: log-compression shells + live tail with
    # mass exp(sum(log_compression)) split over the live walkers
    function ln_evidence_new(df, live, beta)
        w = vcat(ωᵢ(df.log_compression),
                 fill(exp(sum(df.log_compression)) / length(live), length(live)))
        E = vcat(df.emax, live)
        return log(sum(w .* exp.(-beta .* E)))
    end

    # per-N evidence, legacy accounting on the same ledger: iteration-based
    # shells with the conventional ω0 = (K+1)/K, tail q^n split over the live set
    function ln_evidence_legacy(df, live, K, beta)
        n = nrow(df)
        q = K / (K + 1)
        w = vcat(ωᵢ(df.iter, K; ω0=(K + 1) / K),
                 fill(q^n / length(live), length(live)))
        E = vcat(df.emax, live)
        return log(sum(w .* exp.(-beta .* E)))
    end

    # count of bit-exact-zero rows and the longest consecutive run of them
    function zero_tie_stats(df)
        n0 = count(iszero, df.emax)
        best = 0
        cur = 0
        for e in df.emax
            cur = e == 0.0 ? cur + 1 : 0
            best = max(best, cur)
        end
        return n0, best
    end

    # every legal per-cull charge: n/(n+1) for an ordinary cull from n live
    # walkers, (n-1)/n for a plateau eviction; same expressions as src, so
    # membership is bit-exact
    charge_set(K) = Set{Float64}(vcat([log(m / (m + 1)) for m in 1:K],
                                      [log((m - 1) / m) for m in 2:K]))

    ljs_trunc = plateau_cps(rcut)
    ljs_ctrl = plateau_cps(Inf)   # LJParameters accepts cutoff=Inf (its default)
    surface = AtomWalker(surface_system(); freeze_species=[:Ar])
    surface.energy_frozen_part = 0.0u"eV"   # by hand: adsorbate part only
    @test surface.energy_frozen_part == 0.0u"eV"

    E_grid_trunc = field_grid(lj_field(rcut))
    E_grid_ctrl = field_grid(lj_field(Inf))
    f_plateau = plateau_fraction(E_grid_trunc)
    # the untruncated control field has no exact-zero set at all
    @test plateau_fraction(E_grid_ctrl) == 0.0
    # issue #181 measured f = 0.90795 on this same 144^3 grid
    @test isapprox(f_plateau, 0.908; atol=0.002)

    T_gate = 100.0
    beta = 1 / (kB * T_gate)
    ln_z1_trunc = ref_ln_z1(E_grid_trunc, beta)
    ln_z1_ctrl = ref_ln_z1(E_grid_ctrl, beta)

    K = 48
    N1_SEEDS = (11, 12, 13)
    N2_SEEDS = (21, 22, 23)
    CTRL_SEEDS = (41, 42, 43)

    # calibrated gates (see the calibration ledger in the header)
    GATE_EV_N1 = 0.2158
    GATE_MASS_N1 = 0.2573
    GATE_EV_N2 = 0.2647
    GATE_EV_CTRL = 0.3668

    # ---- all RNG use happens here; every gate below is evaluated afterwards ----
    runs_n1 = Dict(s => run_plateau_ns(1, ljs_trunc, surface, s) for s in N1_SEEDS)
    runs_n2 = Dict(s => run_plateau_ns(2, ljs_trunc, surface, s; n_steps=600)
                   for s in N2_SEEDS)
    runs_ctrl = Dict(s => run_plateau_ns(1, ljs_ctrl, surface, s) for s in CTRL_SEEDS)

    # recorded plateau prior mass under both weightings, against the quadrature
    # reference f^N: the accounting statistic the fix changes, whose only
    # sampling noise is the binomial live-set draw of the tie count
    function plateau_mass_devs(df, N)
        zrows = findall(iszero, df.emax)
        shells_new = ωᵢ(df.log_compression)
        P_new = sum(shells_new[zrows])
        shells_leg = ωᵢ(df.iter, K; ω0=(K + 1) / K)
        P_leg = sum(shells_leg[zrows])
        return log(P_new / f_plateau^N), log(P_leg / f_plateau^N), P_leg
    end

    if calib
        for (label, runs, N, lnz) in (("N=1", runs_n1, 1, ln_z1_trunc),
                                      ("N=2", runs_n2, 2, 2 * ln_z1_trunc))
            for s in sort(collect(keys(runs)))
                df, live, po = runs[s]
                dn = ln_evidence_new(df, live, beta) - lnz
                dl = ln_evidence_legacy(df, live, K, beta) - lnz
                dmn, dml, _ = plateau_mass_devs(df, N)
                nz, blk = zero_tie_stats(df)
                println("calib $label seed $s: n_tie=$nz block=$blk ",
                        "live=$(length(live)) refill=$(po.plateau_refill_target)")
                println("  dev_ev_new=$(repr(dn)) dev_ev_leg=$(repr(dl))")
                println("  dev_mass_new=$(repr(dmn)) dev_mass_leg=$(repr(dml))")
            end
        end
        for s in sort(collect(keys(runs_ctrl)))
            df, live, _ = runs_ctrl[s]
            dn = ln_evidence_new(df, live, beta) - ln_z1_ctrl
            println("calib control seed $s: dev_ev_new=$(repr(dn))")
        end
    end

    # ================= item 1: mechanism regression =================
    @testset "mechanism regression: plateau eviction on the issue reproducer" begin
        df, live, params_out = runs_n1[N1_SEEDS[1]]

        # ledger shape: one row per accepted cull, consecutive iters
        @test nrow(df) == 400
        @test df.iter == collect(1:nrow(df))
        @test names(df) == ["iter", "emax", "log_compression"]
        @test all(diff(df.emax) .<= 0)

        # (a) a single block of >= 30 consecutive rows bit-exactly at the
        # plateau energy 0.0 (the issue measures 44-46 at full scale, this
        # reduced scale draws 40-44 across the shipped seeds, expectation
        # K*f = 43.6, binomial sigma = 2.0)
        n_zero, block = zero_tie_stats(df)
        @test block >= 30
        @test block == n_zero          # all zero rows are one consecutive block
        @test block < K                # a whole-live-set tie would fall back
        zrows = findall(iszero, df.emax)
        @test zrows == zrows[1]:zrows[end]
        @test all(df.emax[zrows] .=== 0.0)

        # the eviction schedule charges (n_live - 1)/n_live with the live count
        # shrinking from K, bit-exactly, and the product telescopes to the
        # closed-form block factor (K - n_tie)/K
        @test df.log_compression[zrows] ==
              [log((K - j) / (K - j + 1)) for j in 1:n_zero]
        @test isapprox(exp(sum(df.log_compression[zrows])), (K - n_zero) / K;
                       rtol=1e-12)
        # after the refill the ledger continues below the plateau; the next cull
        # may itself open a secondary micro-tie block from clone energy-residue
        # degeneracy (this seed charges log(47/48) there), so no ordinary-cull
        # assert is made on it beyond legal-charge membership below (guarded: a
        # tie block ending at the final ledger row has no next charge to check)
        if zrows[end] < nrow(df)
            @test df.log_compression[zrows[end]+1] < 0
        end
        # every charge in the ledger is a legal ordinary-cull or eviction factor
        legal = charge_set(K)
        @test all(v -> v in legal, df.log_compression)

        # (b) per-N evidence with the new accounting closes against the
        # midpoint-quadrature reference N*ln z1 at T = 100 K
        dev_ev_new = ln_evidence_new(df, live, beta) - ln_z1_trunc
        @test abs(dev_ev_new) <= GATE_EV_N1
        # and the recorded plateau prior mass closes against ln f
        dev_mass_new, dev_mass_leg, P_leg = plateau_mass_devs(df, 1)
        @test abs(dev_mass_new) <= GATE_MASS_N1

        # (c) the flipped-accounting contrast on the SAME ledger: the legacy
        # iteration-based ωᵢ(iters, K; ω0=(K+1)/K) charges the plateau
        # n_tie*ln((K+1)/K), bounded below 1 nat, instead of ln(1/(1-f)):
        #  - the recorded plateau mass fails the same calibrated gate
        #    (1.87x the gate at calibration, low side);
        @test abs(dev_mass_leg) > GATE_MASS_N1
        @test dev_mass_leg < 0
        #  - equivalently the prior mass assigned below the plateau, which
        #    every evidence and reweighted observable below inherits, is
        #    overestimated: ln(S_leg/S_true) fails high (+1.56 at calibration;
        #    structurally >= +1.37 for any n_tie <= K at this f)
        total_leg = sum(ωᵢ(df.iter, K; ω0=(K + 1) / K)) + (K / (K + 1))^nrow(df)
        S_leg = total_leg - P_leg
        S_true = 1 - f_plateau
        @test log(S_leg / S_true) > 1.0
        #  - and on the same ledger the legacy evidence is strictly biased high
        #    relative to the plateau-aware evidence at any beta > 0 (the plateau
        #    sits at E = 0 and every recorded row below it has E < 0)
        dev_ev_leg = ln_evidence_legacy(df, live, K, beta) - ln_z1_trunc
        @test dev_ev_leg > dev_ev_new

        # (d) the live set is restored to K walkers and the plateau bookkeeping
        # field is cleared
        @test length(live) == K
        @test params_out.plateau_refill_target == 0
    end

    # ================= item 3: exactly solvable per-N closure =================
    @testset "per-N evidence closure, N in (1, 2)" begin
        for (runs, N, seeds, gate_ev, gate_mass, lnz, blk_floor) in
            ((runs_n1, 1, N1_SEEDS, GATE_EV_N1, GATE_MASS_N1, ln_z1_trunc, 30),
             (runs_n2, 2, N2_SEEDS, GATE_EV_N2, nothing, 2 * ln_z1_trunc, 25))
            for s in seeds
                df, live, params_out = runs[s]
                # evidence closure at T = 100 K, per-N calibrated gate
                dev_ev_new = ln_evidence_new(df, live, beta) - lnz
                @test abs(dev_ev_new) <= gate_ev
                # recorded plateau mass vs the quadrature reference f^N;
                # gated at N = 1 only (see the calibration ledger)
                if gate_mass !== nothing
                    dev_mass_new, _, _ = plateau_mass_devs(df, N)
                    @test abs(dev_mass_new) <= gate_mass
                end
                # the same-ledger directional contrast holds for every seed
                dev_ev_leg = ln_evidence_legacy(df, live, K, beta) - lnz
                @test dev_ev_leg > dev_ev_new
                # the initial plateau block is present (expectation K*f^N:
                # 43.6 at N = 1, 39.6 at N = 2)
                n_zero, block = zero_tie_stats(df)
                @test block >= blk_floor
                @test block == n_zero
                # end-state contract: either the run ends outside a tie block
                # with the live set restored, or inside one with the refill
                # target still armed (N = 2 seed 22 ends mid-block)
                @test (params_out.plateau_refill_target == 0 &&
                       length(live) == K) ||
                      (params_out.plateau_refill_target == K &&
                       length(live) < K)
            end
        end
    end

    # ================= scope: untruncated control =================
    @testset "scope: untruncated control (no plateau, uniform tie-free column)" begin
        for s in CTRL_SEEDS
            df, live, params_out = runs_ctrl[s]
            @test nrow(df) == 400
            # zero bit-exact tie rows: no exact zeros, strictly decreasing emax
            n_zero, block = zero_tie_stats(df)
            @test n_zero == 0
            @test all(diff(df.emax) .< 0)
            # the lazily added column is present and uniformly the fixed-K
            # ordinary-cull factor, bit-exactly
            @test names(df) == ["iter", "emax", "log_compression"]
            @test all(df.log_compression .== log(K / (K + 1)))
            # closure at T = 100 K within the control's own calibrated gate
            dev = ln_evidence_new(df, live, beta) - ln_z1_ctrl
            @test abs(dev) <= GATE_EV_CTRL
            # no plateau bookkeeping was ever armed
            @test params_out.plateau_refill_target == 0
            @test length(live) == K
        end
    end

    # ============= scope: lattice ledgers carry no column =============
    @testset "scope: lattice ledgers carry no column (lazy column contract)" begin
        lat = MLattice{1,SquareLattice}(
            lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 1), periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5], components=[[false for _ in 1:16]],
            adsorptions=:full)
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        function fresh_lattice(seed)
            Random.seed!(seed)
            ws = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                  for _ in 1:12]
            for w in ws
                occ = vcat(fill(true, 6), fill(false, 10))
                shuffle!(occ)
                w.configuration.components[1] .= occ
            end
            ls = LatticeGasWalkers(ws, ham; perturb_energy=1e-9)
            p = LatticeNestedSamplingParameters(mc_steps=20,
                energy_perturbation=1e-9, allowed_fail_count=100000)
            return ls, p
        end
        lat_save = SaveEveryN("_test_ties_lat.csv", "_test_ties_lat.traj",
                              "_test_ties_lat.ls", 1000000, 1000000, 1000000)
        ls, p = fresh_lattice(61)
        df, _, _ = nested_sampling(ls, p, 60, MCRandomWalkClone(), lat_save)
        # the lattice step methods keep the four-value return, so their ledgers
        # gain no log_compression column
        @test nrow(df) > 0
        @test names(df) == ["iter", "emax"]
        @test !hasproperty(df, :log_compression)
        # with observables the schema is iter, emax, then the observable columns
        ls, p = fresh_lattice(62)
        df2, _, _ = nested_sampling(ls, p, 60, MCRandomWalkClone(), lat_save;
            observables=[:nocc => (cfg -> Float64(sum(cfg.components[1])))])
        @test nrow(df2) > 0
        @test names(df2) == ["iter", "emax", "nocc"]
        @test !hasproperty(df2, :log_compression)
        rm.(["_test_ties_lat.csv", "_test_ties_lat.traj", "_test_ties_lat.ls"],
           force=true)
    end

    # ============ scope: serial LJAtomWalkers column ============
    @testset "scope: serial LJAtomWalkers column (ledger carries log_compression)" begin
        Random.seed!(51)
        boxp = [[10.0u"Å", 0u"Å", 0u"Å"],
                [0u"Å", 10.0u"Å", 0u"Å"],
                [0u"Å", 0u"Å", 10.0u"Å"]]
        ljp = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0, shift=false)
        mk_walker() = AtomWalker(FastSystem(periodic_system(
            [:Ar => [rand(), rand(), rand()],
             :Ar => [rand(), rand(), rand()],
             :Ar => [rand(), rand(), rand()]], boxp, fractional=true)))
        walkers = [mk_walker() for _ in 1:6]
        ls = LJAtomWalkers(walkers, ljp)
        p = NestedSamplingParameters(mc_steps=100, step_size=0.3,
                                     allowed_fail_count=1000)
        at_save = SaveEveryN(df_filename="_test_ties_at.csv",
                             wk_filename="_test_ties_at.traj.extxyz",
                             ls_filename="_test_ties_at.ls.extxyz",
                             n_traj=10_000_000, n_snap=10_000_000,
                             n_info=10_000_000)
        df, ls_out, params_out = nested_sampling(ls, p, 60, MCRandomWalkClone(),
                                                 at_save)
        rm("_test_ties_at.csv", force=true)
        rm("_test_ties_at.traj.extxyz", force=true)
        rm("_test_ties_at.ls.extxyz", force=true)
        @test nrow(df) > 0
        @test names(df) == ["iter", "emax", "log_compression"]
        @test all(isfinite, df.log_compression)
        @test all(<(0), df.log_compression)
        # uniform-or-tie contract: every charge is a legal ordinary-cull or
        # eviction factor for a live count at most 6
        legal6 = charge_set(6)
        @test all(v -> v in legal6, df.log_compression)
        @test length(ls_out.walkers) == 6
        @test params_out.plateau_refill_target == 0
    end

    # ========== shared fixture: hand-built 6-walker LJAtomWalkers plateau ==========
    # Exercises the OTHER serial step method (plain AtomWalkers; the ledger scope
    # testset above covers it tie-free only) by driving nested_sampling_step!
    # directly through a fully deterministic 4-walker tie block. Cubic all-periodic
    # 30 A box; K = 6 walkers of two H atoms each: 4 walkers hold their pair 12 A
    # apart -- beyond the cutoff radius 3.5*2.5 = 8.75 A, so bit-exactly 0.0 eV --
    # at walker-distinct positions, and 2 walkers hold their pair at the LJ
    # minimum distance (energy < 0). The eviction schedule is then exact:
    # log(5/6), log(4/5), log(3/4), log(2/3), refill on the last step.
    tie6_L = 30.0
    tie6_box = [[tie6_L, 0.0, 0.0]u"Å", [0.0, tie6_L, 0.0]u"Å",
                [0.0, 0.0, tie6_L]u"Å"]
    tie6_lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
    tie6_dmin = 2.0^(1 / 6) * 2.5    # LJ pair-minimum distance, Angstrom

    function tie6_liveset()
        pair(x, y, z, dx, dy) = AtomWalker(FastSystem(atomic_system(
            [:H => [x, y, z]u"Å", :H => [x + dx, y + dy, z]u"Å"],
            tie6_box, (true, true, true))))
        # 4 plateau walkers at distinct positions (distinct z per walker), pair
        # split 12 A along y; 2 sub-plateau walkers at the pair minimum along x
        plateau = [pair(4.0 + 2.0 * i, 6.0, 3.0 + 3.0 * i, 0.0, 12.0)
                   for i in 1:4]
        bound = [pair(15.0, 20.0, 16.0 + 5.0 * i, tie6_dmin, 0.0) for i in 1:2]
        return LJAtomWalkers(vcat(plateau, bound), tie6_lj)
    end

    tie6_params() = NestedSamplingParameters(mc_steps=100, step_size=0.5,
                                             allowed_fail_count=1000)
    TIE6_SCHEDULE = (log(5 / 6), log(4 / 5), log(3 / 4), log(2 / 3))

    all_z(ls) = [ustrip(u"Å", position(w.configuration, i)[3])
                 for w in ls.walkers for i in 1:length(w.configuration)]

    @testset "AtomWalkers serial tie path" begin
        Random.seed!(42)
        ls = tie6_liveset()
        # the constructed plateau is exact: four walkers bit-exactly at 0.0 eV,
        # two strictly below
        @test count(w -> ustrip(u"eV", w.energy) === 0.0, ls.walkers) == 4
        @test count(w -> w.energy < 0.0u"eV", ls.walkers) == 2
        p = tie6_params()
        @test p.plateau_refill_target == 0
        for (k, lt) in enumerate(TIE6_SCHEDULE)
            _, emax, ls, p, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(
                ls, p, MCRandomWalkClone())
            # every step of the block culls at the plateau ceiling, bit-exactly
            @test ustrip(u"eV", emax) === 0.0
            # the eviction schedule charges (n_live - 1)/n_live with the
            # shrinking live count, bit-exactly
            @test log_t === lt
            if k < 4
                # mid-block: the refill target is armed and no refill has run
                @test p.plateau_refill_target == 6
                @test length(ls.walkers) == 6 - k
            end
        end
        # the last tied walker (unique ceiling, target armed) triggered the
        # refill: live set restored, bookkeeping cleared, plateau gone --
        # refilled clones sit strictly below it
        @test length(ls.walkers) == 6
        @test p.plateau_refill_target == 0
        @test all(w -> w.energy < 0.0u"eV", ls.walkers)
    end

    @testset "2D-dims refill stays on the constrained manifold" begin
        Random.seed!(43)
        ls = tie6_liveset()
        @test count(w -> ustrip(u"eV", w.energy) === 0.0, ls.walkers) == 4
        @test count(w -> w.energy < 0.0u"eV", ls.walkers) == 2
        # the z multiset of the live set before any step; 2D walks and 2D
        # refills never move z (2D dims are only supported on the plain
        # AtomWalkers method, exactly the method under test)
        z0 = Set(all_z(ls))
        p = tie6_params()
        for lt in TIE6_SCHEDULE
            _, emax, ls, p, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(
                ls, p, MCRandomWalkClone(dims=[1, 2]))
            @test ustrip(u"eV", emax) === 0.0
            @test log_t === lt
        end
        @test length(ls.walkers) == 6
        @test p.plateau_refill_target == 0
        @test all(w -> w.energy < 0.0u"eV", ls.walkers)
        # every z in the refilled live set is exactly a member of the original
        # z multiset: survivors keep theirs, 2D-refilled clones inherit their
        # ancestor's (a 3D refill regression would move z off the manifold)
        @test all(z -> z in z0, all_z(ls))
    end
end
