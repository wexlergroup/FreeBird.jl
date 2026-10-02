@testset "Lattice walks that carry their tie-breaking offsets" begin
    using Random

    # The lattice walks break ties between arrangements of equal energy with a
    # perturbation drawn uniformly on [-delta/2, delta/2). By default every
    # proposal draws a new one (`perturbation_mode = :per_proposal`), so inside
    # a level of equal energy a proposal passes the ceiling only if its new
    # perturbation lands below the ceiling's: d units of depth into the level,
    # with probability about exp(-d). Under `perturbation_mode = :carried` a
    # walker tests each proposal with its own perturbation and redraws it after
    # every step from its law below the ceiling, so moves that keep the
    # ordering variable pass at any depth. The model below is the
    # one-dimensional lattice gas with nearest-neighbor attraction on a
    # periodic chain of 12 sites (the Ising chain in lattice-gas form).

    cp_M = 12
    cp_eps = -0.02
    cp_ring() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(cp_M, 1, 1),
        periodicity=(true, false, false), cutoff_radii=[1.1],
        components=[fill(false, cp_M)], adsorptions=:full)
    cp_ham = GenericLatticeHamiltonian(0.0, [cp_eps], u"eV")
    cp_ham_meV = GenericLatticeHamiltonian(0.0, [1000 * cp_eps], u"meV")
    cp_raw(w, h) = ustrip(u"eV", interacting_energy(w.configuration, h))
    cp_n(w) = count(w.configuration.components[1])

    # A walker holding the given occupied sites, with its energy set to the
    # unperturbed energy plus `offset` (eV)
    function cp_walker(sites, h; offset=0.0)
        lattice = cp_ring()
        lattice.components[1][sites] .= true
        w = LatticeWalker(lattice)
        w.energy = (cp_raw(w, h) + offset) * u"eV"
        return w
    end
    # A ceiling `depth` units of depth into the level at unperturbed key `key`:
    # delta*exp(-depth) above the bottom of the offset range, whose part below
    # the ceiling then holds exp(-depth) of the range
    cp_ceiling(key, delta, depth) = key - delta / 2 + delta * exp(-depth)

    @testset "modes and constructors" begin
        @test MCMixedMoves().perturbation_mode === :per_proposal
        @test MCMixedMoves(5, 1).perturbation_mode === :per_proposal
        @test MCMixedMoves(1, 0, 0, 0.3, 0.3, 50, 0.01, 1.0).perturbation_mode === :per_proposal
        @test MCMixedMoves(1, 0, 0, 0.3, 0.3, 50, 0.01, 1.0, true).perturbation_mode === :per_proposal
        loose = MCMixedMoves(1, 0, 0, 0.3, 0.3, 50, 0.01, 1, true)   # an Int for a Float64 field converts
        @test loose.cluster_p_ceiling === 1.0 && loose.perturbation_mode === :per_proposal
        @test MCMixedMoves(perturbation_mode=:carried).perturbation_mode === :carried
        @test MCGrandCanonicalMoves().perturbation_mode === :per_proposal
        @test MCGrandCanonicalMoves(perturbation_mode=:carried).perturbation_mode === :carried
        @test_throws ArgumentError MCMixedMoves(perturbation_mode=:fresh)
        @test_throws ArgumentError MCGrandCanonicalMoves(perturbation_mode=:fresh)
        w = cp_walker([1, 4, 7], cp_ham)
        @test_throws ArgumentError MC_random_walk!(1, w, cp_ham, 1.0; perturbation_mode=:fresh)
        @test_throws ArgumentError MC_cluster_walk!(1, w, cp_ham, 1.0, 0.3; perturbation_mode=:fresh)
        @test_throws ArgumentError MC_grand_canonical_walk!(1, w, cp_ham, 1.0, 0.0;
                                                            perturbation_mode=:fresh)
    end

    @testset "the refreshed offset follows its law below the ceiling" begin
        refresh! = FreeBird.MonteCarloMoves._refresh_carried_offset!
        delta = 1e-3
        for h in (cp_ham, cp_ham_meV), (mu, sites) in ((0.0, [1, 2, 7]), (-0.03, [1, 2, 7]))
            w = cp_walker(sites, h)
            raw = interacting_energy(w.configuration, h)          # the Hamiltonian's unit
            key = cp_raw(w, h) - mu * cp_n(w)
            ceiling = cp_ceiling(key, delta, 3.0)                # inside the level
            top = ceiling - key                                  # the offsets below the ceiling end here
            offsets = Float64[]
            for u in (0.0, 0.25, 0.5, 0.75, 0.999)
                new = refresh!(w, raw, 0.0u"eV", u, ceiling * u"eV", delta, mu, cp_n(w))
                push!(offsets, ustrip(u"eV", new))
                @test w.energy.val - mu * cp_n(w) < ceiling
                @test w.energy.val ≈ cp_raw(w, h) + offsets[end] atol = 1e-15
            end
            @test offsets[1] ≈ -delta / 2 atol = 1e-15
            @test issorted(offsets)
            @test offsets[3] ≈ (-delta / 2 + top) / 2 atol = 1e-12
            @test offsets[end] < top
            # at the top of the range the new key can round onto the ceiling: then the
            # offset is kept, so the key never reaches the ceiling
            new = refresh!(w, raw, offsets[end] * u"eV", prevfloat(1.0), ceiling * u"eV", delta, mu, cp_n(w))
            @test ustrip(u"eV", new) < top
            @test w.energy.val - mu * cp_n(w) < ceiling
            @test w.energy.val ≈ cp_raw(w, h) + ustrip(u"eV", new) atol = 1e-15
            # an arrangement well below the ceiling draws from the whole range
            new = refresh!(w, raw, 0.0u"eV", prevfloat(1.0), (key + 1.0) * u"eV", delta, mu, cp_n(w))
            @test ustrip(u"eV", new) ≈ delta / 2 atol = 1e-12
        end
        # no offset below the ceiling, or delta = 0: the offset and the walker are kept
        w = cp_walker([1, 2, 7], cp_ham; offset=-4e-4)
        raw = interacting_energy(w.configuration, cp_ham)
        e = w.energy
        @test refresh!(w, raw, -4e-4u"eV", 0.5, (cp_raw(w, cp_ham) - 1e-3) * u"eV", 1e-3, 0.0, 3) == -4e-4u"eV"
        @test w.energy == e
        @test refresh!(w, raw, -4e-4u"eV", 0.5, (cp_raw(w, cp_ham) + 1.0) * u"eV", 0.0, 0.0, 3) == -4e-4u"eV"
        @test w.energy == e
        # a negative delta spans the same range, as the default mode's perturbation does
        w1, w2 = cp_walker([1, 2, 7], cp_ham), cp_walker([1, 2, 7], cp_ham)
        c = (cp_raw(w1, cp_ham) + 1.0) * u"eV"
        @test refresh!(w1, raw, 0.0u"eV", 0.3, c, -1e-3, 0.0, 3) == refresh!(w2, raw, 0.0u"eV", 0.3, c, 1e-3, 0.0, 3)
        # with a meV Hamiltonian, redraws near the top of narrow ranges never store a key at
        # or above the ceiling, the key computed in the walker's unit as the drivers do
        stored = above = 0
        for sites in ([1, 2, 7], [1, 2, 3, 9], [2, 3, 5, 6, 10]), mu in (0.0, -0.013, -0.03), d in (1e-6, 1e-12)
            w = cp_walker(sites, cp_ham_meV)
            raw = interacting_energy(w.configuration, cp_ham_meV)
            n = length(sites)
            key = cp_raw(w, cp_ham_meV) - mu * n
            for k in 1:400, u in (prevfloat(1.0), 1 - 1e-9, 1 - 1e-12, 0.999999999988668)
                c = key - d / 2 + d * k / 401 * 1e-3
                w.energy = (cp_raw(w, cp_ham_meV) - d / 2) * u"eV"
                before = w.energy
                refresh!(w, raw, (-d / 2) * u"eV", u, c * u"eV", d, mu, n)
                if w.energy != before
                    stored += 1
                    above += !(w.energy.val - mu * n < c)
                end
            end
        end
        @test stored > 1000
        @test above == 0
    end

    @testset "moves that keep the ordering variable at depth 20 into its level" begin
        delta = 1e-3
        depth = 20.0          # a new offset lands below the ceiling with probability exp(-20) = 2e-9

        # Fixed N: three particles in one block, the lowest level for N = 3. The
        # proposals that keep it are the swaps of two sites of equal occupancy
        # (90 of the 144 ordered site pairs) and the exchanges of an end particle
        # with the empty site past the other end (4 of 144: either order of a pair
        # swaps the same two sites), so a walk that is not frozen accepts about
        # 94 of 144, 65%
        block = [1, 2, 3]
        key = cp_raw(cp_walker(block, cp_ham), cp_ham)
        ceiling = cp_ceiling(key, delta, depth)
        for (mode, incremental) in ((:per_proposal, false), (:per_proposal, true),
                                    (:carried, false), (:carried, true))
            w = cp_walker(block, cp_ham; offset=-delta / 2 + delta * exp(-depth) / 2)
            Random.seed!(2026)
            accepted, rate, w = MC_random_walk!(2000, w, cp_ham, ceiling; energy_perturb=delta,
                                                incremental=incremental, perturbation_mode=mode)
            if mode === :per_proposal
                @test !accepted
                @test findall(w.configuration.components[1]) == block
            else
                @test rate > 0.5
                @test cp_raw(w, cp_ham) ≈ key atol = 1e-12  # still one block of three
                @test w.energy.val < ceiling
            end
        end

        # Grand-canonical walk ordered by U - mu*N at mu = eps: every single block has
        # U - mu*N = -eps, so insertions and deletions at a block's ends keep the key
        mu = cp_eps
        sites = collect(1:6)
        key = cp_raw(cp_walker(sites, cp_ham), cp_ham) - mu * length(sites)
        ceiling = cp_ceiling(key, delta, depth)
        for mode in (:per_proposal, :carried)
            w = cp_walker(sites, cp_ham; offset=-delta / 2 + delta * exp(-depth) / 2)
            Random.seed!(2027)
            accepted, rate, w, _, _, stats = MC_grand_canonical_walk!(2000, w, cp_ham, ceiling, mu;
                p_move=0.2, p_insert=0.4, energy_perturb=delta, perturbation_mode=mode)
            changed = stats.insert_uniform_accepted + stats.delete_accepted
            if mode === :per_proposal
                @test !accepted
                @test changed == 0
                @test findall(w.configuration.components[1]) == sites
            else
                # a block that shrinks to nothing or fills the chain drops below the ceiling's
                # level and stays there, which takes at least six changes of its length
                @test changed >= 6
                @test w.energy.val - mu * cp_n(w) < ceiling
            end
        end

        # Reference-fugacity walk (mu = 0, z0 = 0.01) from the empty chain, the ceiling
        # inside U = 0: an insertion that forms no pair keeps U = 0
        ceiling = cp_ceiling(0.0, delta, depth)
        for mode in (:per_proposal, :carried)
            w = cp_walker(Int[], cp_ham; offset=-delta / 2 + delta * exp(-depth) / 2)
            Random.seed!(2028)
            accepted, rate, w, _, _, stats = MC_grand_canonical_walk!(2000, w, cp_ham, ceiling, 0.0;
                p_move=0.4, p_insert=0.3, z0=0.01, energy_perturb=delta, perturbation_mode=mode)
            if mode === :per_proposal
                @test !accepted
                @test stats.insert_uniform_accepted == 0
                @test cp_n(w) == 0
            else
                @test stats.insert_uniform_accepted > 0      # about 50 on average, 14 at the least in 1000 seeds
                @test w.energy.val < ceiling
            end
        end
    end

    @testset "the lattice nested-sampling steps pass the mode to the walks" begin
        # Walkers in one level of equal ordering key with their offsets packed 20 units of depth into it: the clone's walk
        # accepts a move that keeps the key only under :carried, and fail_count records whether it moved
        delta, depth, K = 1e-3, 20.0, 8
        packed(i) = -delta / 2 + delta * exp(-depth) * i / (K + 1)
        live(arrs) = LatticeGasWalkers([cp_walker(arrs[mod1(i, length(arrs))], cp_ham; offset=packed(i)) for i in 1:K],
                                       cp_ham; assign_energy=false)
        for mode in (:per_proposal, :carried)
            moved = mode === :carried ? 0 : 1
            for (wf, cf) in ((1, 0), (0, 1))    # swaps alone, clusters alone (fixed N = 3, one block)
                params = LatticeNestedSamplingParameters(mc_steps=50, energy_perturbation=delta, allowed_fail_count=10^9)
                Random.seed!(2030)
                nested_sampling_step!(live([[1, 2, 3], [5, 6, 7], [9, 10, 11]]), params,
                                      MCMixedMoves(walks_freq=wf, clusters_freq=cf, perturbation_mode=mode))
                @test params.fail_count == moved
            end
            # ordered by U - mu*N at mu = eps: single blocks of four to seven sites
            params = GrandCanonicalNestedSamplingParameters(mc_steps=50, chemical_potential=cp_eps,
                                                            energy_perturbation=delta, allowed_fail_count=10^9)
            Random.seed!(2031)
            nested_sampling_step!(live([collect(1:4), collect(3:8), collect(5:11)]), params,
                                  MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, perturbation_mode=mode))
            @test params.fail_count == moved
            # reference fugacity from the empty chain (U = 0)
            params = IdealGasReferencedGCNSParameters(mc_steps=50, reference_fugacity=0.5, energy_perturbation=delta,
                                                      allowed_fail_count=10^9)
            Random.seed!(2032)
            nested_sampling_step!(live([Int[]]), params,
                                  MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, perturbation_mode=mode))
            @test params.fail_count == moved
        end
    end

    @testset "after every walk the stored energy is the arrangement's energy plus a redrawn offset below the ceiling" begin
        # A deterministic invariant: from a random arrangement with the ceiling inside its own level, short walks of every
        # kind; afterwards the stored key lies strictly below the ceiling, the stored energy minus the arrangement's energy
        # (the offset) lies in [-delta/2, min(delta/2, ceiling - key)), and the offset was redrawn (the start is neither
        # empty nor full, so the first step is never a guard skip or a biased null proposal, which redraw nothing)
        delta = 1e-3
        Random.seed!(2050)
        bad_key = bad_offset = not_redrawn = 0
        walks = 0
        for trial in 1:200
            sites = Int[]
            while !(0 < length(sites) < cp_M)
                sites = findall(rand(cp_M) .< rand())
            end
            mu = rand((0.0, cp_eps, -0.05))
            w0 = cp_walker(sites, cp_ham)
            key0 = cp_raw(w0, cp_ham) - mu * cp_n(w0)
            c = key0 - delta / 2 + delta * rand()                 # inside the level
            off0 = -delta / 2 + rand() * (c - key0 + delta / 2)
            for kind in 1:6
                mu != 0.0 && kind <= 3 && continue                # the fixed-N walks are energy-ordered
                w = cp_walker(sites, cp_ham; offset=off0)
                n_steps = rand(1:3)
                if kind == 1
                    MC_random_walk!(n_steps, w, cp_ham, c; energy_perturb=delta, perturbation_mode=:carried)
                elseif kind == 2
                    MC_random_walk!(n_steps, w, cp_ham, c; energy_perturb=delta, incremental=true, perturbation_mode=:carried)
                elseif kind == 3
                    MC_cluster_walk!(n_steps, w, cp_ham, c, 0.5; energy_perturb=delta, perturbation_mode=:carried)
                else
                    kw = kind == 4 ? (;) : kind == 5 ? (incremental=true,) : (clusters_freq=1, p_bias=0.5)
                    MC_grand_canonical_walk!(n_steps, w, cp_ham, c, mu; p_move=0.4, p_insert=0.3, z0=rand((1.0, 0.2)),
                                             energy_perturb=delta, perturbation_mode=:carried, kw...)
                end
                key = cp_raw(w, cp_ham) - mu * cp_n(w)
                off = w.energy.val - cp_raw(w, cp_ham)
                walks += 1
                bad_key += !(w.energy.val - mu * cp_n(w) < c)
                bad_offset += !(-delta / 2 - 1e-12 <= off < min(delta / 2, c - key) + 1e-12)
                not_redrawn += off == off0
            end
        end
        @test walks > 500
        @test bad_key == 0
        @test bad_offset == 0
        @test not_redrawn == 0
    end

    @testset "the carried walks leave the law restricted below the ceiling invariant" begin
        # Walkers drawn exactly from the prior restricted below a ceiling that lies inside a level (a quarter of the level's
        # offset range below it), walked 1 and 50 steps each: the fraction in the ceiling's level must stay at its exact
        # value (a 5-sigma binomial band; one step catches a wrong first test, 50 steps a wrong redraw schedule)
        delta = 1e-3
        function exact_draws(arrs, keys, level, n)
            c = level - delta / 2 + 0.25 * delta
            width = [clamp(min(delta / 2, c - k) + delta / 2, 0.0, delta) for k in keys]
            p = sum(width[abs.(keys .- level) .< 1e-9]) / sum(width)
            cum = cumsum(width) ./ sum(width)
            draws = map(1:n) do _
                j = min(searchsortedfirst(cum, rand()), length(arrs))
                (arrs[j], -delta / 2 + rand() * width[j])
            end
            return c, p, draws
        end
        masks = 0:(2^cp_M - 1)
        sites_of(m) = findall(digits(m, base=2, pad=cp_M) .== 1)
        n = 8000
        # fixed N = 3 (one block, U = -2 eps; pair and single, U = -eps; apart, U = 0), the ceiling in U = eps
        arrs = [sites_of(m) for m in masks if count_ones(m) == 3]
        keys = [cp_raw(cp_walker(a, cp_ham), cp_ham) for a in arrs]
        for kind in 1:3, steps in (1, 50)
            Random.seed!(2060 + kind + steps)
            c, p, draws = exact_draws(arrs, keys, cp_eps, n)
            inlevel = 0
            for (a, off) in draws
                w = cp_walker(a, cp_ham; offset=off)
                if kind == 3
                    MC_cluster_walk!(steps, w, cp_ham, c, 0.5; energy_perturb=delta, perturbation_mode=:carried)
                else
                    MC_random_walk!(steps, w, cp_ham, c; energy_perturb=delta, incremental=kind == 2,
                                    perturbation_mode=:carried)
                end
                inlevel += abs(cp_raw(w, cp_ham) - cp_eps) < 1e-9
            end
            @test abs(inlevel / n - p) < 5 * sqrt(p * (1 - p) / n)
        end
        # ordered by U - mu*N at mu = eps (single blocks at -eps, the empty and the full chain at 0), z0 = 1
        arrs = [sites_of(m) for m in masks]
        keys = [cp_raw(cp_walker(a, cp_ham), cp_ham) - cp_eps * length(a) for a in arrs]
        for (i, kw) in enumerate(((;), (incremental=true,), (clusters_freq=1, p_bias=0.5))), steps in (1, 50)
            Random.seed!(2070 + i + steps)
            c, p, draws = exact_draws(arrs, keys, -cp_eps, n)
            inlevel = 0
            for (a, off) in draws
                w = cp_walker(a, cp_ham; offset=off)
                MC_grand_canonical_walk!(steps, w, cp_ham, c, cp_eps; p_move=0.4, p_insert=0.3, energy_perturb=delta,
                                         perturbation_mode=:carried, kw...)
                inlevel += abs(cp_raw(w, cp_ham) - cp_eps * cp_n(w) + cp_eps) < 1e-9
            end
            @test abs(inlevel / n - p) < 5 * sqrt(p * (1 - p) / n)
        end
    end

    @testset "every step consumes the same number of draws in either mode" begin
        # One-step walks from identical walkers and seeds: the generator's state after the step does not depend on the mode
        delta = 1e-3
        Random.seed!(2080)
        states = [(findall(rand(cp_M) .< f), s) for f in (0.0, 0.3, 0.6, 1.0) for s in (-0.4, 0.2, 2.0)]
        mismatches = 0
        for (sites, s) in states, kind in 1:5, seed in 1:20
            after = map((:per_proposal, :carried)) do mode
                w = cp_walker(sites, cp_ham; offset=-delta / 4)
                c = w.energy.val - (kind >= 4 ? cp_eps * cp_n(w) : 0.0) + s * delta
                Random.seed!(seed)
                if kind == 1
                    MC_random_walk!(1, w, cp_ham, c; energy_perturb=delta, perturbation_mode=mode)
                elseif kind == 2
                    MC_random_walk!(1, w, cp_ham, c; energy_perturb=delta, incremental=true, perturbation_mode=mode)
                elseif kind == 3
                    MC_cluster_walk!(1, w, cp_ham, c, 0.5; energy_perturb=delta, perturbation_mode=mode)
                else
                    MC_grand_canonical_walk!(1, w, cp_ham, c, cp_eps; p_move=0.4, p_insert=0.3, energy_perturb=delta,
                                             perturbation_mode=mode, (kind == 5 ? (clusters_freq=1, p_bias=0.5) : (;))...)
                end
                copy(Random.default_rng())
            end
            mismatches += after[1] != after[2]
        end
        @test mismatches == 0
    end

    @testset "agreement with exact enumeration under :carried" begin
        # The three lattice nested-sampling drivers with perturbation_mode = :carried
        # on a 16-site chain, against the exact enumeration of its 2^16 arrangements
        # at 300 K: the mean over seeds 1-3 of each run's error in ln Xi (ln Z_N for
        # fixed N). A run's estimate counts its discarded and final walkers with the
        # geometric shells of the default compression (the i-th discarded walker
        # stands for (e^(1/K) - 1) e^(-i/K) of the prior, each final walker for
        # e^(-n/K)/K). Tolerances: three times the largest |mean| of ten disjoint
        # seed triples (seeds 1-30) measured on this code: 0.109 for the reference
        # fugacity, 0.214 for fixed N, 0.131 for fixed mu. An integration check of
        # the drivers under :carried: at these settings it does not separate the
        # modes (over 600 seeds both read about +0.01 for each driver), which the
        # step tests below do.
        M, K, walk, delta, T, kB = 16, 20, 64, 1e-6, 300.0, 8.617333262e-5
        chain() = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(M, 1, 1),
            periodicity=(true, false, false), cutoff_radii=[1.1],
            components=[fill(false, M)], adsorptions=:full)
        lse(x) = (t = maximum(x); t + log(sum(exp.(x .- t))))
        logC(n, k) = sum(log, (n - k + 1):n; init=0.0) - sum(log, 1:k; init=0.0)
        shells(n) = vcat(log(expm1(1 / K)) .- (1:n) ./ K, fill(-n / K - log(K), K))
        no_files = SaveEveryN(df_filename="unused.csv", wk_filename="unused.extxyz",
                              ls_filename="unused.ls.extxyz", n_traj=typemax(Int),
                              n_snap=typemax(Int), n_info=typemax(Int))
        exact_N, exact_U = Int[], Float64[]
        let lattice = chain()
            for mask in 0:(2^M - 1)
                for site in 1:M
                    lattice.components[1][site] = ((mask >> (site - 1)) & 1) == 1
                end
                push!(exact_N, count(lattice.components[1]))
                push!(exact_U, ustrip(u"eV", interacting_energy(lattice, cp_ham)))
            end
        end
        exact_lnXi(mu) = lse((mu .* exact_N .- exact_U) ./ (kB * T))
        live_N(live) = [count(w.configuration.components[1]) for w in live.walkers]
        live_U(live) = [w.energy.val for w in live.walkers]

        # reference fugacity z0 = e^-2, at its own chemical potential kT ln z0
        errors = map(1:3) do seed
            z0 = exp(-2.0)
            walkers = LatticeGasWalkers(replicate_walkers(chain(), K), cp_ham; assign_energy=false)
            params = IdealGasReferencedGCNSParameters(mc_steps=walk, reference_fugacity=z0,
                energy_perturbation=delta, allowed_fail_count=10^9)
            Random.seed!(seed)
            steps = round(Int, (1.15 * M * log1p(1 / z0) + 2) * K)
            df, live, _ = ideal_gas_referenced_nested_sampling(walkers, params, steps,
                MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, perturbation_mode=:carried), no_files)
            N, U = vcat(df.num_particles, live_N(live)), vcat(df.emax, live_U(live))
            mu = -2.0 * kB * T
            lse(M * log1p(z0) .+ shells(size(df, 1)) .+ 2.0 .* N .+ (mu .* N .- U) ./ (kB * T)) - exact_lnXi(mu)
        end
        @test abs(sum(errors) / 3) < 0.35

        # fixed N = 8
        errors = map(1:3) do seed
            Random.seed!(seed)
            walkers = replicate_walkers(chain(), K)
            for w in walkers
                w.configuration.components[1][randperm(M)[1:8]] .= true
            end
            live = LatticeGasWalkers(walkers, cp_ham; assign_energy=true, perturb_energy=delta)
            params = LatticeNestedSamplingParameters(mc_steps=walk, energy_perturbation=delta,
                                                     allowed_fail_count=10^9)
            steps = round(Int, (1.15 * logC(M, 8) + 2) * K)
            df, live, _ = nested_sampling(live, params, steps,
                MCMixedMoves(walks_freq=1, clusters_freq=0, perturbation_mode=:carried), no_files)
            U = vcat(df.emax, live_U(live))
            lse(logC(M, 8) .+ shells(size(df, 1)) .- U ./ (kB * T)) - lse(-exact_U[exact_N .== 8] ./ (kB * T))
        end
        @test abs(sum(errors) / 3) < 0.65

        # fixed mu = eps, where the empty and the full chain are the two lowest arrangements
        errors = map(1:3) do seed
            walkers = LatticeGasWalkers(replicate_walkers(chain(), K), cp_ham; assign_energy=false)
            params = GrandCanonicalNestedSamplingParameters(mc_steps=walk, chemical_potential=cp_eps,
                energy_perturbation=delta, allowed_fail_count=10^9)
            Random.seed!(seed)
            steps = round(Int, (1.15 * M * log(2) + 2) * K)
            df, live, _ = grand_canonical_nested_sampling(walkers, params, steps,
                MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, perturbation_mode=:carried), no_files)
            N, U = vcat(df.num_particles, live_N(live)), vcat(df.energy, live_U(live))
            lse(M * log(2) .+ shells(size(df, 1)) .+ (cp_eps .* N .- U) ./ (kB * T)) - exact_lnXi(cp_eps)
        end
        @test abs(sum(errors) / 3) < 0.40
    end
end
