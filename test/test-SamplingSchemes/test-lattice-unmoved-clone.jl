@testset "Lattice NS steps keep an unmoved clone" begin
    using Random

    kb = 8.617333262e-5  # eV/K

    # A replacement walker is a copy of a surviving walker, moved by a walk
    # that stays below the ceiling. When the walk accepts no move, the copy
    # still holds its parent's arrangement, which lies below the ceiling, so
    # the lattice steps keep it and redraw only its tie-breaking offset from
    # the offset's law below the ceiling (`_keep_unmoved_clone!`). Retrying
    # from another parent instead would copy walkers that cannot move too
    # rarely. `mc_steps = 0` below makes every walk accept no move.

    # Attractive nearest-neighbor lattice gas on a periodic 3 x 3 square lattice
    uc_lattice() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(3, 3, 1),
        periodicity=(true, true, false), cutoff_radii=[1.1],
        components=[fill(false, 9)], adsorptions=:full)
    uc_ham = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
    uc_raw(w) = ustrip(u"eV", interacting_energy(w.configuration, uc_ham))
    uc_n(w) = count(w.configuration.components[1])
    uc_moves = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3, incremental=true,
                                     swap_mode=:occupied_empty)

    # K walkers with seeded random occupations (exactly N particles when given)
    function uc_liveset(K, delta; N=nothing, seed=1)
        Random.seed!(seed)
        ws = replicate_walkers(uc_lattice(), K)
        for w in ws
            occ = w.configuration.components[1]
            occ .= false
            occ[randperm(9)[1:(N === nothing ? rand(0:9) : N)]] .= true
        end
        return LatticeGasWalkers(ws, uc_ham; assign_energy=true,
                                 perturb_energy=delta)
    end

    @testset "ideal-gas-referenced step" begin
        delta = 1e-3
        ls = uc_liveset(6, delta; seed=91)
        sort_by_energy!(ls)
        ceiling = ls.walkers[1].energy.val
        worst_n = uc_n(ls.walkers[1])
        others = [copy(w.configuration.components[1]) for w in ls.walkers[2:end]]
        params = IdealGasReferencedGCNSParameters(mc_steps=0,
            reference_fugacity=0.1, energy_perturbation=delta)
        iter, emax, n_worst, ls2, params2 = nested_sampling_step!(ls, params, uc_moves)
        @test iter isa Int
        @test emax.val == ceiling
        @test n_worst == worst_n
        @test length(ls2.walkers) == 6
        kept = ls2.walkers[end]
        @test kept.configuration.components[1] in others
        @test kept.energy.val < ceiling
        @test abs(kept.energy.val - uc_raw(kept)) <= delta / 2
        @test params2.fail_count == 1
    end

    @testset "grand-canonical (Omega-sorted) step" begin
        delta = 1e-3
        mu = -0.05
        ls = uc_liveset(6, delta; seed=92)
        omega(w) = w.energy.val - mu * uc_n(w)
        ceiling = maximum(omega(w) for w in ls.walkers)
        worst = argmax(omega, ls.walkers)
        others = [copy(w.configuration.components[1]) for w in ls.walkers if w !== worst]
        params = GrandCanonicalNestedSamplingParameters(mc_steps=0,
            chemical_potential=mu, energy_perturbation=delta)
        iter, omega_max, energy, n_worst, ls2, params2 =
            nested_sampling_step!(ls, params, uc_moves)
        @test iter isa Int
        @test omega_max.val == ceiling
        @test n_worst == uc_n(worst)
        @test length(ls2.walkers) == 6
        kept = ls2.walkers[end]
        @test kept.configuration.components[1] in others
        @test omega(kept) < ceiling
        @test abs(kept.energy.val - uc_raw(kept)) <= delta / 2
        @test params2.fail_count == 1
    end

    @testset "fixed-N steps ($(nameof(typeof(routine))))" for routine in (
            MCMixedMoves(walks_freq=1, clusters_freq=0, incremental=true),
            MCRandomWalkClone())
        delta = 1e-3
        ls = uc_liveset(6, delta; N=4, seed=93)
        sort_by_energy!(ls)
        ceiling = ls.walkers[1].energy.val
        others = [copy(w.configuration.components[1]) for w in ls.walkers[2:end]]
        params = LatticeNestedSamplingParameters(mc_steps=0, energy_perturbation=delta)
        iter, emax, ls2, params2 = nested_sampling_step!(ls, params, routine)
        @test iter isa Int
        @test emax.val == ceiling
        @test length(ls2.walkers) == 6
        kept = ls2.walkers[end]
        @test kept.configuration.components[1] in others
        @test kept.energy.val < ceiling
        @test abs(kept.energy.val - uc_raw(kept)) <= delta / 2
        @test params2.fail_count == 1
    end

    @testset "MCRandomWalkMaxE still records nothing" begin
        # This routine walks a copy of the culled walker itself, which is not
        # below the ceiling, so an unmoved copy is not kept
        delta = 1e-3
        ls = uc_liveset(6, delta; N=4, seed=94)
        before = sort([w.energy.val for w in ls.walkers])
        params = LatticeNestedSamplingParameters(mc_steps=0, energy_perturbation=delta)
        iter, emax, ls2, params2 = nested_sampling_step!(ls, params, MCRandomWalkMaxE())
        @test iter === missing
        @test ismissing(emax)
        @test sort([w.energy.val for w in ls2.walkers]) == before
        @test params2.fail_count == 1
    end

    @testset "an arrangement tied with the ceiling is not kept" begin
        # Without tie-breaking offsets and with every walker in one arrangement,
        # no replacement lies strictly below the ceiling: the step records
        # nothing, as before
        Random.seed!(95)
        ws = replicate_walkers(uc_lattice(), 5)
        for w in ws
            w.configuration.components[1][[1, 2, 5]] .= true
        end
        ls = LatticeGasWalkers(ws, uc_ham; assign_energy=true, perturb_energy=0.0)
        params = GrandCanonicalNestedSamplingParameters(mc_steps=0,
            chemical_potential=-0.05, energy_perturbation=0.0)
        out = nested_sampling_step!(ls, params, uc_moves)
        @test out[1] === missing
        @test length(out[5].walkers) == 5
        @test params.fail_count == 1
    end

    @testset "a stall stop keeps the row of the step that triggered it ($(name))" for (name, run) in (
            ("ideal-gas-referenced", (ls, n) -> ideal_gas_referenced_nested_sampling(ls,
                IdealGasReferencedGCNSParameters(mc_steps=0, reference_fugacity=0.1,
                    energy_perturbation=1e-3, allowed_fail_count=3), n, uc_moves,
                SaveEveryN("test_uc.csv", "test_uc.traj", "test_uc.ls", 100000, 100000, 100000);
                stop_on_stall=true)),
            ("Omega-sorted", (ls, n) -> grand_canonical_nested_sampling(ls,
                GrandCanonicalNestedSamplingParameters(mc_steps=0, chemical_potential=-0.05,
                    energy_perturbation=1e-3, allowed_fail_count=3), n, uc_moves,
                SaveEveryN("test_uc.csv", "test_uc.traj", "test_uc.ls", 100000, 100000, 100000);
                stop_on_stall=true)))
        # Every walk accepts no move (mc_steps = 0) and keeps its clone, so
        # each step culls and records a row while fail_count counts up; the
        # third step reaches allowed_fail_count and stops the run after its
        # row is recorded
        ls = LatticeGasWalkers(replicate_walkers(uc_lattice(), 6), uc_ham; assign_energy=false)
        Random.seed!(98)
        df, ls2, params = run(ls, 20)
        rm.(["test_uc.csv", "test_uc.traj", "test_uc.ls"], force=true)
        @test nrow(df) == 3
        @test df.iter == [1, 2, 3]
        @test length(ls2.walkers) == 6
        @test params.fail_count == 3
    end

    @testset "the redrawn offset follows its law below the ceiling" begin
        keep! = FreeBird.SamplingSchemes._keep_unmoved_clone!
        delta = 1e-3
        w = uc_liveset(1, 0.0; N=3, seed=96).walkers[1]
        u = uc_raw(w)
        Random.seed!(97)
        # ceiling inside the offset range: offsets uniform on [-delta/2, 0.2 delta)
        offsets = Float64[]
        for _ in 1:4000
            @assert keep!(w, uc_ham, u + 0.2 * delta, delta)
            push!(offsets, w.energy.val - u)
        end
        @test all(-delta / 2 .<= offsets .< 0.2 * delta)
        @test abs(mean(offsets) - (-0.15 * delta)) < 5 * 0.7 * delta / sqrt(12 * 4000)
        # ceiling above the range: the whole range
        offsets = [(keep!(w, uc_ham, u + delta, delta); w.energy.val - u) for _ in 1:4000]
        @test all(-delta / 2 .<= offsets .< delta / 2)
        @test abs(mean(offsets)) < 5 * delta / sqrt(12 * 4000)
        # no offset below the ceiling: not kept, energy unchanged
        e_before = w.energy
        @test !keep!(w, uc_ham, u - delta / 2, delta)
        @test w.energy == e_before
        # the U - mu N key
        mu = -0.05
        n = uc_n(w)
        @test keep!(w, uc_ham, u - mu * n + 0.1 * delta, delta; mu=mu)
        @test w.energy.val - mu * n < u - mu * n + 0.1 * delta
        # without offsets: kept, at the unperturbed energy, only strictly below the ceiling
        @test keep!(w, uc_ham, u + 1e-12, 0.0)
        @test w.energy.val == u
        @test !keep!(w, uc_ham, u, 0.0)
        # random numbers: one draw when the offset range reaches below the
        # ceiling, none without offsets or when the range lies above it
        Random.seed!(7)
        r = rand(UInt64)
        Random.seed!(7)
        keep!(w, uc_ham, u + 1e-12, 0.0)
        @test rand(UInt64) == r
        Random.seed!(7)
        keep!(w, uc_ham, u - delta, delta)
        @test rand(UInt64) == r
        Random.seed!(7)
        rand()
        r = rand(UInt64)
        Random.seed!(7)
        keep!(w, uc_ham, u + delta, delta)
        @test rand(UInt64) == r
    end

    @testset "a kept clone lies strictly below the ceiling in floating point" begin
        # With the ceiling a few units in the last place above the lowest
        # offset, U + offset can round onto the ceiling; the walks reject a
        # proposal at the ceiling, and such a clone is not kept either
        keep! = FreeBird.SamplingSchemes._keep_unmoved_clone!
        w = uc_liveset(1, 0.0; N=3, seed=96).walkers[1]
        n = uc_n(w)
        Random.seed!(101)
        ok = Bool[]
        kept = 0
        for mu in (0.0, -0.05), delta in (1e-12, 1e-6), k in (1, 4)
            ceiling = (uc_raw(w) - mu * n) - delta / 2
            for _ in 1:k
                ceiling = nextfloat(ceiling)
            end
            for _ in 1:1000
                before = w.energy
                if keep!(w, uc_ham, ceiling, delta; mu=mu)
                    kept += 1
                    push!(ok, w.energy.val - mu * n < ceiling)
                else
                    push!(ok, w.energy == before)
                end
            end
        end
        @test all(ok)
        @test kept > 0
    end

    @testset "energies in the walker's unit (meV Hamiltonians)" begin
        # The lattice walks compare energies with the ceiling, and store them,
        # in the walker's unit (eV) whatever the Hamiltonian's unit; so does
        # the kept clone, with attractive and repulsive couplings alike
        keep! = FreeBird.SamplingSchemes._keep_unmoved_clone!
        delta = 1e-3
        w = uc_liveset(1, 0.0; N=3, seed=96).walkers[1]
        Random.seed!(102)
        for h in (GenericLatticeHamiltonian(-40.0, [-10.0], u"meV"),
                  GenericLatticeHamiltonian(40.0, [10.0], u"meV"))
            u = ustrip(u"eV", interacting_energy(w.configuration, h))
            offsets = [(keep!(w, h, u + 0.2 * delta, delta) ? w.energy.val - u : NaN)
                       for _ in 1:2000]
            @test all(-delta / 2 .<= offsets .< 0.2 * delta)
            @test abs(mean(offsets) - (-0.15 * delta)) < 5 * 0.7 * delta / sqrt(12 * 2000)
        end
    end

    @testset "ideal-gas-referenced sampler at small z0 against exact enumeration" begin
        # At ln z0 = -8 the empty lattice holds 99.7% of the prior and a walk
        # from it moves only by an insertion accepted with probability
        # 9 z0 per attempt, so most walks from it accept no move. Calibration
        # (seeds 1 to 12, this fixture): with the unmoved clones kept, the
        # largest mean of any three seeds deviates from the exact values by
        # 0.058 in <N> and in ln Xi (gates at 0.2); when every walk that
        # accepts no move is instead retried from another parent, seeds 1 to 3
        # deviate by +0.68 in <N> (exact 0.069) and +0.55 in ln Xi.
        K, walk, lnz0, T = 50, 64, -8.0, 150.0
        mu = kb * T * lnz0
        steps = round(Int, (1.15 * 9 * log1p(exp(-lnz0)) + 2) * K)
        # exact reference: all 2^9 arrangements
        lat = uc_lattice()
        occ = lat.components[1]
        logw = Float64[]
        ns = Int[]
        for code in 0:(2^9 - 1)
            for s in 1:9
                occ[s] = isodd(code >> (s - 1))
            end
            push!(ns, count(occ))
            push!(logw, (mu * count(occ) - ustrip(u"eV", interacting_energy(lat, uc_ham))) / (kb * T))
        end
        top = maximum(logw)
        exact_lnXi = top + log(sum(exp.(logw .- top)))
        exact_N = sum(exp.(logw .- exact_lnXi) .* ns)
        save = SaveEveryN("test_uc.csv", "test_uc.traj", "test_uc.ls", 100000, 100000, 100000)
        dN, dX = Float64[], Float64[]
        for seed in 1:3
            ls = LatticeGasWalkers(replicate_walkers(uc_lattice(), K), uc_ham; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=walk,
                reference_fugacity=exp(lnz0), energy_perturbation=1e-6,
                allowed_fail_count=10^9)
            Random.seed!(seed)
            df, lsx, _ = ideal_gas_referenced_nested_sampling(ls, p, steps, uc_moves, save)
            rm.(["test_uc.csv", "test_uc.traj", "test_uc.ls"], force=true)
            stats = gc_thermodynamic_stats_ideal_ref(df, 9, exp(lnz0), [mu], [T], K;
                live_emax=[w.energy.val for w in lsx.walkers],
                live_numbers=[count(w.configuration.components[1]) for w in lsx.walkers])
            push!(dN, stats.mean_N[1, 1] - exact_N)
            push!(dX, stats.logXi[1, 1] - exact_lnXi)
        end
        @test abs(mean(dN)) < 0.2
        @test abs(mean(dX)) < 0.2
    end
end
