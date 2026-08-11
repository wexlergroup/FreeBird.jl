@testset "Hard-core lattice models (finite-J athermal recipe)" begin
    using Random

    # Hard squares (square lattice, NN exclusion) and hard hexagons
    # (triangular lattice, NN exclusion) via the finite-J recipe: a finite
    # repulsive coupling J on a single-shell lattice makes the energy
    # J × (number of excluded-neighbor pairs), an integer-leveled ladder whose
    # E = 0 manifold is exactly the hard-core configuration space. Post-
    # processing at β·J ≥ (ladder depth in nats) + ~40 turns the Boltzmann
    # factor into an exact 0/1 indicator, so every observable is athermal —
    # a function of the activity z = exp(βμ) alone. See the hard-core recipe
    # in the GenericLatticeHamiltonian docstring.

    kb = 8.617333262e-5  # eV/K
    J = 1.0              # eV per excluded-neighbor pair
    ham_hc = GenericLatticeHamiltonian(0.0, [J], u"eV")

    hc_square() = MLattice{1,SquareLattice}(
        lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)],
        supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false),
        cutoff_radii=[1.1],
        components=:equal,
        adsorptions=:full)

    # (3, 3, 1) with the 2-site basis: M = 18, commensurate with the
    # √3×√3 three-sublattice ordering (close packing N = M/3 = 6, g = 3).
    hc_tri() = MLattice{1,TriangularLattice}(
        lattice_constant=1.0,
        supercell_dimensions=(3, 3, 1),
        periodicity=(true, true, false),
        cutoff_radii=[1.1],
        components=:equal,
        adsorptions=:full)

    function lattice_at(make, N::Int)
        lat = make()
        lat.components[1] .= false
        lat.components[1][1:N] .= true
        return lat
    end

    # Brute-force reference built from the neighbor lists alone (independent
    # of EnergyEval): bitmask adjacency, count independent sets by N.
    function neighbor_masks(lat)
        M = length(lat.components[1])
        masks = zeros(UInt64, M)
        for i in 1:M
            m = UInt64(0)
            for j in lat.neighbors[i][1]
                m |= UInt64(1) << (j - 1)
            end
            masks[i] = m
        end
        return masks
    end

    function independent_set_counts(masks, M)
        g = zeros(Int, M + 1)
        for cfg in UInt64(0):(UInt64(1) << M) - UInt64(1)
            ok = true
            c = cfg
            while c != 0
                i = trailing_zeros(c) + 1
                if masks[i] & cfg != 0
                    ok = false
                    break
                end
                c &= c - UInt64(1)
            end
            ok && (g[count_ones(cfg) + 1] += 1)
        end
        return g
    end

    function violation_count(masks, occ)
        v = 0
        for i in eachindex(occ)
            occ[i] || continue
            for j in eachindex(occ)
                if occ[j] && masks[i] & (UInt64(1) << (j - 1)) != 0
                    v += 1
                end
            end
        end
        return v ÷ 2
    end

    function triangle_count(masks, M)
        t = 0
        for i in 1:M, j in (i + 1):M
            masks[i] & (UInt64(1) << (j - 1)) == 0 && continue
            common = masks[i] & masks[j] & ~((UInt64(1) << j) - UInt64(1))
            t += count_ones(common)
        end
        return t
    end

    # Exact hard-core grand-canonical stats from the g_N table:
    # Ξ(z) = Σ_N g_N z^N.
    function exact_hc_stats(g, z)
        Ns = [N for N in 0:length(g)-1 if g[N+1] > 0]
        lt = [log(g[N+1]) + N * log(z) for N in Ns]
        mx = maximum(lt)
        w = exp.(lt .- mx)
        sw = sum(w)
        return (logXi=mx + log(sw),
                mean_N=sum(w .* Ns) / sw,
                var_N=sum(w .* Ns .^ 2) / sw - (sum(w .* Ns) / sw)^2)
    end

    sq = hc_square()
    tr = hc_tri()
    M_sq = length(sq.components[1])
    M_tr = length(tr.components[1])
    masks_sq = neighbor_masks(sq)
    masks_tr = neighbor_masks(tr)
    g_sq = independent_set_counts(masks_sq, M_sq)
    g_tr = independent_set_counts(masks_tr, M_tr)

    save_strategy = SaveEveryN(
        df_filename="test_hc_df.csv",
        wk_filename="test_hc.traj.extxyz",
        ls_filename="test_hc.ls.extxyz",
        n_traj=1_000_000, n_snap=1_000_000, n_info=1_000_000)
    hc_cleanup() = foreach(f -> rm(f, force=true),
        ("test_hc_df.csv", "test_hc.traj.extxyz", "test_hc.ls.extxyz"))

    # ================================================================
    @testset "exact enumeration references" begin
        # Hard squares, 4×4 PBC: pins the square neighbor lists and the
        # close-packed c(2×2) pair.
        @test M_sq == 16
        @test all(length(sq.neighbors[i][1]) == 4 for i in 1:M_sq)
        @test g_sq[1:9] == [1, 16, 88, 208, 228, 128, 56, 16, 2]
        @test all(g_sq[10:end] .== 0)

        # Hard hexagons, 18-site commensurate torus: pins the triangular
        # neighbor lists, the cell commensurability, and the three √3×√3
        # sublattice states at close packing.
        @test M_tr == 18
        @test all(length(tr.neighbors[i][1]) == 6 for i in 1:M_tr)
        @test g_tr[1:7] == [1, 18, 99, 180, 99, 18, 3]
        @test all(g_tr[8:end] .== 0)

        # Pin the adjacency graphs exactly. The square torus is bipartite
        # (no triangles). The 18-site cell is a genuine quotient of the
        # triangular lattice, but its x-circumference is 3, so besides the
        # 2M = 36 face triangles it carries 6 wrap-around 3-cycles (one per
        # a₁-row): 42 total. The exact reference is enumerated on this same
        # graph, so the pipeline comparison is exact for this cell; it is a
        # pipeline validation, not a thermodynamic-limit statement —
        # production cells must use circumference ≥ 4.
        @test triangle_count(masks_sq, M_sq) == 0
        @test triangle_count(masks_tr, M_tr) == 42
    end

    # ================================================================
    @testset "energy = J × violation count" begin
        Random.seed!(7)
        for (make, masks, M) in ((hc_square, masks_sq, M_sq),
                                 (hc_tri, masks_tr, M_tr))
            lat = make()
            for _ in 1:50
                occ = rand(Bool, M)
                lat.components[1] .= occ
                @test interacting_energy(lat, ham_hc) ≈
                      J * violation_count(masks, occ) * u"eV"
            end
        end
    end

    # ================================================================
    # Fixed-N route: per-N canonical NS ladders descend the violation count
    # to the E = 0 manifold; gc_thermodynamic_stats_fixed_N assembles Ξ(z).
    # This is the primary route for athermal models: every z is an exact
    # polynomial evaluation in the per-N evidence (no reweighting window).
    function fixed_N_hard_core(make, N_max; K, n_iters, mc_steps, seed,
                               routine=MCRandomWalkClone())
        Random.seed!(seed)
        dfs = Vector{DataFrame}(undef, N_max + 1)
        live = Vector{Vector{Float64}}(undef, N_max + 1)
        dfs[1] = DataFrame(iter=Int[], emax=Float64[])
        live[1] = Float64[]
        for N in 1:N_max
            walkers = [
                begin
                    lat = lattice_at(make, N)
                    generate_random_new_lattice_sample!(lat)
                    LatticeWalker(lat)
                end for _ in 1:K]
            liveset = LatticeGasWalkers(walkers, ham_hc, perturb_energy=1e-9)
            ns_params = LatticeNestedSamplingParameters(
                mc_steps=mc_steps,
                energy_perturbation=1e-9,
                allowed_fail_count=100_000)
            df, final_ls, _ = nested_sampling(
                liveset, ns_params, n_iters, routine, save_strategy)
            dfs[N + 1] = df
            live[N + 1] = [ustrip(u"eV", w.energy) for w in final_ls.walkers]
        end
        hc_cleanup()
        return dfs, live
    end

    @testset "fixed-N route vs exact" begin
        K = 60
        # Ladder depth 900/60 = 15 nats; the athermal evaluation temperature
        # must satisfy β·J ≥ depth + ~40: T = 150 K gives β·J ≈ 77.
        n_iters = 900
        T_ath = 150.0
        zs = [0.5, 1.0, 2.0, 3.796, 11.09]  # numeric z_c of hard squares and
                                            # Baxter's z_c = (11+5√5)/2 among them
        μ_grid = [kb * T_ath * log(z) for z in zs] .* u"eV"

        for (tag, make, M, N_max, g, seed) in (
                ("hard squares", hc_square, M_sq, 8, g_sq, 1000),
                ("hard hexagons", hc_tri, M_tr, 6, g_tr, 2000))
            @testset "$tag" begin
                dfs, live = fixed_N_hard_core(make, N_max;
                    K=K, n_iters=n_iters, mc_steps=100, seed=seed)

                # Every ladder descended into the allowed (E = 0) manifold:
                # the athermal termination gate.
                for N in 1:N_max
                    @test issorted(dfs[N + 1].emax, rev=true)
                    @test all(e -> e < J / 2, live[N + 1])
                end

                stats = gc_thermodynamic_stats_fixed_N(
                    dfs, collect(0:N_max), M, μ_grid, [T_ath * u"K"];
                    n_walkers=K, n_cull=1, ω0=(K + 1) / K, live_emax=live)

                # Tolerances sized ≥ 3σ of the intrinsic NS shrinkage noise:
                # the deepest sectors carry H_N = ln(C(M,N)/g_N) ≈ 8.8 nats,
                # so σ(log Z_N) ≈ √(H_N/K) ≈ 0.38 per sector and
                # σ(logΞ) ≈ 0.23 at the largest z after the sector mixture.
                # atol = 0.75 ≈ 3.3σ; still discriminating (a dropped
                # (1+z0)^M-style normalization shifts logΞ by ≫ 1).
                for (k, z) in enumerate(zs)
                    ex = exact_hc_stats(g, z)
                    @test isapprox(stats.logXi[k, 1], ex.logXi, atol=0.75)
                    @test isapprox(stats.mean_N[k, 1], ex.mean_N, atol=0.5)
                    # Athermal: ⟨E⟩ is pure tie-breaking noise
                    @test abs(stats.mean_U[k, 1]) < 1e-6
                end

                # Per-N evidence: at the athermal temperature log_Z_N is the
                # log independent-set count log g_N. atol = 1.5 ≈ 3.9σ of the
                # per-sector noise, ≪ the log C(M,N) ≥ 7.5 shift a dropped
                # binomial factor would cause.
                @test stats.N_values == collect(0:N_max)
                @test size(stats.log_Z_N) == (N_max + 1, 1)
                for N in 0:N_max
                    @test isapprox(stats.log_Z_N[N + 1, 1], log(g[N + 1]),
                                   atol=1.5)
                end
            end
        end
    end

    # ================================================================
    # Cluster-move A/B: the triangular geometric cluster move is variance
    # reduction only, so mixed-move ladders must reproduce the same exact
    # per-N evidence as the local-move reference above. Unbiasedness is the
    # assertion; a hard variance comparison over a few seeds would be
    # F-statistics noise and is deliberately not asserted.
    @testset "fixed-N route with cluster moves (hard hexagons)" begin
        K = 60
        n_iters = 900
        T_ath = 150.0
        routine = MCMixedMoves(walks_freq=1, clusters_freq=1)
        dfs, live = fixed_N_hard_core(hc_tri, 6;
            K=K, n_iters=n_iters, mc_steps=100, seed=4000, routine=routine)

        # Every mixed-move ladder still descends into the E = 0 manifold
        for N in 1:6
            @test issorted(dfs[N + 1].emax, rev=true)
            @test all(e -> e < J / 2, live[N + 1])
        end

        stats = gc_thermodynamic_stats_fixed_N(
            dfs, collect(0:6), M_tr, [0.0] .* u"eV", [T_ath * u"K"];
            n_walkers=K, n_cull=1, ω0=(K + 1) / K, live_emax=live)

        # Same per-sector tolerance as the local-move reference (≥ 3σ)
        for N in 0:6
            @test isapprox(stats.log_Z_N[N + 1, 1], log(g_tr[N + 1]), atol=1.5)
        end

        # Wiring: the mixed step exercised cluster moves, tuned cluster_p
        # adaptively, and conserved N through the full NS machinery
        Random.seed!(4001)
        walkers = [
            begin
                lat = lattice_at(hc_tri, 4)
                generate_random_new_lattice_sample!(lat)
                LatticeWalker(lat)
            end for _ in 1:20]
        liveset = LatticeGasWalkers(walkers, ham_hc, perturb_energy=1e-9)
        ns_params = LatticeNestedSamplingParameters(
            mc_steps=40, energy_perturbation=1e-9, allowed_fail_count=100_000)
        _, ls_out, params_out = nested_sampling(
            liveset, ns_params, 200, routine, save_strategy)
        hc_cleanup()
        @test !isempty(params_out.cluster_p_history)
        @test all(sum(w.configuration.components[1]) == 4 for w in ls_out.walkers)
    end

    # ================================================================
    # igref route: one grand-canonical run against the z0-Bernoulli prior;
    # the descent through the violating shells measures the allowed set's
    # prior mass, and the (1+z0)^M normalization stays valid because the
    # prior support is never restricted.
    @testset "igref route vs exact (hard squares, z0 = 1)" begin
        Random.seed!(7)
        K = 100
        template = hc_square()
        template.components[1] .= false
        walkers = [LatticeWalker(deepcopy(template), energy=0.0u"eV", iter=0)
                   for _ in 1:K]
        liveset = LatticeGasWalkers(walkers, ham_hc; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=1.0,
            energy_perturbation=1e-9, allowed_fail_count=100_000)
        mc_routine = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)
        df, final_ls, _ = ideal_gas_referenced_nested_sampling(
            liveset, params, Int64(1500), mc_routine, save_strategy)
        hc_cleanup()

        @test nrow(df) > 0
        @test issorted(df.emax, rev=true)
        # The ladder reached the allowed manifold and the live set sits in it
        @test minimum(df.emax) < J / 2
        live_E = [w.energy.val for w in final_ls.walkers]
        live_N = [sum(w.configuration.components[1]) for w in final_ls.walkers]
        @test all(e -> e < J / 2, live_E)

        T_ath = 150.0
        zs = [0.7, 1.0, 1.5, 2.0]
        μs = [kb * T_ath * log(z) for z in zs]
        stats = gc_thermodynamic_stats_ideal_ref(
            df, M_sq, 1.0, μs, [T_ath], K;
            ω0=(K + 1) / K, live_emax=live_E, live_numbers=live_N)

        # Tolerances ≥ 3σ of the NS evidence noise (σ(logΞ) ≈ 0.25 at this
        # depth/K), while a dropped (1+z0)^M normalization would shift logΞ
        # by M·ln 2 ≈ 11 — still cleanly discriminating.
        for (i, z) in enumerate(zs[1:3])  # z = 2.0 is the N_eff check below
            ex = exact_hc_stats(g_sq, z)
            @test isapprox(stats.logXi[i, 1], ex.logXi, atol=0.75)
            @test isapprox(stats.mean_N[i, 1], ex.mean_N, atol=0.5)
        end

        # Kish N_eff degrades away from the run's z0 — the trust-window
        # diagnostic every reweighted grid point must be gated on.
        @test stats.N_eff[4, 1] < stats.N_eff[2, 1]
    end

    # ================================================================
    @testset "igref route vs exact (hard hexagons, z0 = 1)" begin
        # The 18-site cell's close-packed states are exactly the three pure
        # √3×√3 sublattice states: zero energy and maximal order parameter.
        # A deterministic tie between g_tr's close-packing degeneracy 3 and
        # the sublattice geometry, independent of the sampled run below.
        let lat = hc_tri()
            for c in 0:2
                lat.components[1] .= [
                    begin
                        b = (s - 1) % 2
                        ci = ((s - 1) ÷ 2) % 3
                        (b == 0 ? ci % 3 : (ci + 2) % 3) == c
                    end for s in 1:M_tr]
                @test interacting_energy(lat, ham_hc) ≈ 0.0u"eV" atol=1e-12u"eV"
                @test order_parameter_sqrt3(lat) == 1 / 3
            end
        end

        # Full-compression depth M·ln(1 + z0) = 18·ln 2 ≈ 12.5 nats at K = 100
        # wants ≈ 1.15·K·D ≈ 1440 iterations; 1800 adds margin (the square
        # testset's 1500 at D ≈ 11.1, scaled).
        Random.seed!(22)
        K = 100
        template = hc_tri()
        template.components[1] .= false
        walkers = [LatticeWalker(deepcopy(template), energy=0.0u"eV", iter=0)
                   for _ in 1:K]
        liveset = LatticeGasWalkers(walkers, ham_hc; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=1.0,
            energy_perturbation=1e-9, allowed_fail_count=100_000)
        mc_routine = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)
        df, final_ls, _ = ideal_gas_referenced_nested_sampling(
            liveset, params, Int64(1800), mc_routine, save_strategy)
        hc_cleanup()

        @test nrow(df) > 0
        @test issorted(df.emax, rev=true)
        # The ladder reached the allowed manifold and the live set sits in it
        @test minimum(df.emax) < J / 2
        live_E = [w.energy.val for w in final_ls.walkers]
        live_N = [sum(w.configuration.components[1]) for w in final_ls.walkers]
        @test all(e -> e < J / 2, live_E)

        # Any close-packed live walker must be one of the three sublattice
        # states (the only N = 6 members of the allowed manifold). At z0 = 1
        # their prior weight is 3/418 per walker, so this is usually
        # vacuously true; the deterministic check above is the real tie.
        for w in final_ls.walkers
            sum(w.configuration.components[1]) == 6 || continue
            @test order_parameter_sqrt3(w.configuration) == 1 / 3
        end

        T_ath = 150.0
        zs = [0.7, 1.0, 1.5, 2.0]
        μs = [kb * T_ath * log(z) for z in zs]
        stats = gc_thermodynamic_stats_ideal_ref(
            df, M_tr, 1.0, μs, [T_ath], K;
            ω0=(K + 1) / K, live_emax=live_E, live_numbers=live_N)

        # Tolerances ≥ 3σ of the evidence noise, sized from three-seed
        # scatter at this configuration (max observed |ΔlogΞ| 0.37,
        # |Δ⟨N⟩| 0.24); a dropped (1+z0)^M normalization would shift logΞ
        # by 18·ln 2 ≈ 12.5 — still cleanly discriminating.
        for (i, z) in enumerate(zs[1:3])  # z = 2.0 is the N_eff check below
            ex = exact_hc_stats(g_tr, z)
            @test isapprox(stats.logXi[i, 1], ex.logXi, atol=0.75)
            @test isapprox(stats.mean_N[i, 1], ex.mean_N, atol=0.5)
        end

        # Kish N_eff degrades away from the run's z0
        @test stats.N_eff[4, 1] < stats.N_eff[2, 1]
    end

    # ================================================================
    @testset "biased :cavity walk samples the restricted prior (hard squares)" begin
        # Fixed ceiling J/2 in (0, J) with mu = 0, z0 = 1: the only reachable
        # states are independent sets (E = J x violations), and the walk's
        # stationary law is the RESTRICTED Bernoulli prior
        # P(N) = g_N / sum(g) = g_N / 743.
        # A fixed-ceiling walk is used instead of an NS dead-point histogram:
        # dead points are moving-ceiling order statistics with no closed-form
        # N-marginal, while this chain has an exact target, and it exercises
        # the composite kernel directly, including the isolated-particle
        # delete branch (on an independent set every occupied site is
        # isolated, so every deletion takes that branch).
        Random.seed!(158)
        template_bc = hc_square()
        template_bc.components[1] .= false
        walker_bc = LatticeWalker(template_bc, energy=0.0u"eV", iter=0)
        wcap = J / 2
        gc_walk!(n) = MC_grand_canonical_walk!(n, walker_bc, ham_hc, wcap, 0.0;
            p_move=0.4, p_insert=0.3, z0=1.0,
            p_bias=0.9, bias_predicate=:cavity, bias_shells=1)
        gc_walk!(500)   # burn-in from the empty (allowed) state

        n_samples = 12_000
        counts = zeros(Int, M_sq + 1)
        rate_sum = 0.0
        all_independent = true
        for _ in 1:n_samples
            _, r, _, _, _, _ = gc_walk!(15)   # 15-step decorrelation blocks
            occ = walker_bc.configuration.components[1]
            all_independent &= violation_count(masks_sq, occ) == 0
            counts[sum(occ) + 1] += 1
            rate_sum += r
        end

        # The ceiling was never breached: every sample is an independent set
        @test all_independent
        @test all(counts[10:end] .== 0)   # g_N = 0 for N > 8

        # Multinomial gates: p_N = g_N/743, sigma = sqrt(p(1-p)/n) at
        # n = 12000, with mild residual autocorrelation from the 15-step
        # blocks. Three-seed calibration (seeds 158/159/160): max |dev| =
        # 2.3 iid-sigma, zero tail counts, acceptance rate 0.70 on every
        # seed; the shipped k = 7 is >= 3x the maximum. Discrimination: an
        # uncorrected biased kernel distorts P(N) at the tens-of-percent
        # level, orders above these gates.
        Zg = sum(g_sq)   # = 743
        for N in 0:8
            p = g_sq[N + 1] / Zg
            sig = sqrt(p * (1 - p) / n_samples)
            @test abs(counts[N + 1] / n_samples - p) < 7 * sig
        end

        # Acceptance sanity: cavity insertions never create a violation, so
        # the ceiling never rejects the biased branch; only the MH factor
        # can. A soft floor pins that the kernel is not ceiling-starved.
        @test rate_sum / n_samples > 0.1
    end
end
