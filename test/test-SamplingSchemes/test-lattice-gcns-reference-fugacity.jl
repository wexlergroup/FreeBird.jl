@testset "Lattice GC-NS reference fugacity" begin
    using Random
    using Unitful
    using DataFrames

    # The reference_fugacity keyword of GrandCanonicalNestedSamplingParameters
    # selects the ideal-lattice-gas prior of the Omega-sorted lattice driver
    # (per-site fugacity z0: a configuration with N particles carries prior
    # weight z0^N, total mass (1 + z0)^M, the uniform prior at z0 = 1), draws
    # the initial live set from it, and passes z0 to the insert/delete
    # acceptances of MC_grand_canonical_walk!. The matching keyword of
    # gc_thermodynamic_stats restores the flat counting measure with a
    # per-shell factor z0^(-N_j). Every run here is seeded; the one
    # statistical gate (the exact-enumeration testset) reuses the calibrated
    # precedent of test-grand_canonical_ns.jl and says so where it applies.

    rf_lat1() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false), cutoff_radii=[1.1],
        components=[[false for _ in 1:16]], adsorptions=:full)
    rf_ham1() = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
    rf_lat2() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false), cutoff_radii=[1.1, 1.5],
        components=[[false for _ in 1:16]], adsorptions=:full)
    rf_ham2() = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

    rf_save(tag) = SaveEveryN("t_rf_$(tag).csv", "t_rf_$(tag).traj",
                              "t_rf_$(tag).ls", 1000000, 1000000, 1000000)
    rf_cleanup(tag) = rm.(["t_rf_$(tag).csv", "t_rf_$(tag).traj",
                           "t_rf_$(tag).ls"], force=true)

    rf_keys = (:swap_attempted, :swap_accepted,
        :cluster_attempted, :cluster_accepted,
        :insert_uniform_attempted, :insert_uniform_accepted,
        :insert_biased_attempted, :insert_biased_accepted,
        :delete_attempted, :delete_accepted)

    rf_live_E(ls) = [w.energy.val for w in ls.walkers]
    rf_live_N(ls) = [sum(w.configuration.components[1]) for w in ls.walkers]

    @testset "constructor" begin
        gc = GrandCanonicalNestedSamplingParameters()
        @test gc.reference_fugacity == 1.0
        @test gc.init_occupation_p === nothing
        gc2 = GrandCanonicalNestedSamplingParameters(reference_fugacity=0.5)
        @test gc2.reference_fugacity == 0.5
        @test gc2.init_occupation_p === nothing
        # the positivity check of the energy-sorted sibling
        @test_throws ArgumentError GrandCanonicalNestedSamplingParameters(
            reference_fugacity=0.0)
        @test_throws ArgumentError GrandCanonicalNestedSamplingParameters(
            reference_fugacity=-1.0)
        # a number for init_occupation_p is compatible with the uniform prior
        # only: the reference measure fixes the initial law otherwise
        @test_throws ArgumentError GrandCanonicalNestedSamplingParameters(
            init_occupation_p=0.3, reference_fugacity=0.5)
        gc3 = GrandCanonicalNestedSamplingParameters(init_occupation_p=0.3)
        @test gc3.init_occupation_p == 0.3
        @test gc3.reference_fugacity == 1.0
        gc4 = GrandCanonicalNestedSamplingParameters(init_occupation_p=0.3,
                                                     reference_fugacity=1.0)
        @test gc4.init_occupation_p == 0.3
    end

    @testset "initial live set: prior draw and counter reset" begin
        # 2000 walkers x 16 sites at z0 = 0.5: the per-site occupancy is a
        # Bernoulli(1/3) sample of size 32000, sigma 0.0026; the gate atol
        # 0.02 is about 7.6 sigma. Stale iteration counters are zeroed.
        Random.seed!(90341)
        ws = replicate_walkers(rf_lat1(), 2000)
        for w in ws
            w.iter = 7
        end
        ls = LatticeGasWalkers(ws, rf_ham1(); assign_energy=false)
        gc = GrandCanonicalNestedSamplingParameters(reference_fugacity=0.5,
                                                    energy_perturbation=1e-9)
        SamplingSchemes._init_gc_walkers!(ls, gc)
        @test isapprox(sum(rf_live_N(ls)) / (2000 * 16), 1 / 3; atol=0.02)
        @test all(w.iter == 0 for w in ls.walkers)
        # a number for init_occupation_p under the uniform prior keeps the
        # previous free choice of the initial law
        Random.seed!(90342)
        ls3 = LatticeGasWalkers(replicate_walkers(rf_lat1(), 2000), rf_ham1();
                                assign_energy=false)
        SamplingSchemes._init_gc_walkers!(ls3,
            GrandCanonicalNestedSamplingParameters(init_occupation_p=0.3,
                                                   energy_perturbation=1e-9))
        @test isapprox(sum(rf_live_N(ls3)) / (2000 * 16), 0.3; atol=0.02)
    end

    @testset "default identity (same-process replay)" begin
        # The keyword's default, an explicit 1.0, and the previous default
        # spelled out (init_occupation_p = 0.5, the value z0/(1.0 + z0) takes
        # at z0 = 1.0 exactly) replay the same seeded trajectory: the kernel's
        # z0 default was 1.0 already, and the initial draw is the same
        # rand() < 0.5 per site. The absolute pins in
        # test-lattice-gc-trajectory-pins.jl guard the same path against the
        # pre-keyword trajectory.
        function rf_run_omega(seed, tag; kwargs...)
            ls = LatticeGasWalkers(replicate_walkers(rf_lat1(), 10), rf_ham1();
                                   assign_energy=false)
            gc = GrandCanonicalNestedSamplingParameters(; mc_steps=30,
                chemical_potential=-0.05, energy_perturbation=1e-9, kwargs...)
            Random.seed!(seed)
            df, lsx, pout = grand_canonical_nested_sampling(ls, gc, Int64(50),
                MCGrandCanonicalMoves(), rf_save(tag))
            rf_cleanup(tag)
            return df, lsx, pout
        end
        dA, lsA, pA = rf_run_omega(90301, "a")
        dB, lsB, pB = rf_run_omega(90301, "b"; reference_fugacity=1.0)
        dC, lsC, pC = rf_run_omega(90301, "c"; init_occupation_p=0.5)
        @test nrow(dA) > 0
        @test dA.iter == dB.iter == dC.iter
        @test dA.omega == dB.omega == dC.omega
        @test dA.energy == dB.energy == dC.energy
        @test dA.num_particles == dB.num_particles == dC.num_particles
        @test rf_live_E(lsA) == rf_live_E(lsB) == rf_live_E(lsC)
        for k in rf_keys
            @test pA.move_stats[k] == pB.move_stats[k] == pC.move_stats[k]
        end
        # The init redraws every occupancy, reassigns the energy and zeroes
        # the counter, so a live set that already ran (counters at the
        # previous run's depth) replays the fresh-walker ledger under the
        # same seed; without the counter reset the iter column would carry
        # the previous run's offset
        Random.seed!(90302)
        dR, _, _ = grand_canonical_nested_sampling(lsA,
            GrandCanonicalNestedSamplingParameters(mc_steps=30,
                chemical_potential=-0.05, energy_perturbation=1e-9),
            Int64(50), MCGrandCanonicalMoves(), rf_save("r"))
        rf_cleanup("r")
        dF, _, _ = rf_run_omega(90302, "f")
        @test dR.iter == dF.iter
        @test dR.omega == dF.omega
        @test dR.energy == dF.energy
        @test dR.num_particles == dF.num_particles
    end

    @testset "mu = 0 identity with the energy-sorted driver" begin
        # At chemical_potential = 0 the Omega-sorted driver and the
        # energy-sorted driver (ideal_gas_referenced_nested_sampling) run the
        # same algorithm on the same prior, so the same seed must replay the
        # same ledger. The identity rests on four conditions:
        #   1. fresh walkers built the same way for each call (energy 0, iter
        #      0): both inits redraw every occupancy with rand() < z0/(1 + z0),
        #      assign the perturbed energy with one rand(), and zero the
        #      counter, in that order, so they consume the stream identically;
        #   2. energy_perturbation != 0 (the energy-sorted constructor refuses
        #      zero; distinct perturbed keys also make the tie order moot);
        #   3. Random.seed! immediately before each driver call;
        #   4. the same MCGrandCanonicalMoves instance, mc_steps, n_steps and
        #      save frequencies.
        # Per iteration both steps then sort the same keys descending (at
        # mu = 0 the cached Omega key E - 0.0 * N is the energy bit for bit,
        # and the cached sortperm reproduces sort!(by = energy) including the
        # tie order, as pinned in test-grand_canonical_ns.jl), draw the parent
        # with one rand over the same eligible set, and call
        # MC_grand_canonical_walk! with the same ceiling (omega_worst.val is
        # emax.val), mu = 0.0, z0 and n_max (the Omega-sorted default
        # typemax(Int64) is the kernel's own default), so every draw
        # coincides. The Omega-sorted omega column then equals its energy
        # column, which equals the energy-sorted emax column.
        z0 = 0.3
        moves = MCGrandCanonicalMoves()
        rf_fresh() = LatticeGasWalkers(replicate_walkers(rf_lat1(), 12),
                                       rf_ham1(); assign_energy=false)
        lsO = rf_fresh()
        gc = GrandCanonicalNestedSamplingParameters(mc_steps=30,
            chemical_potential=0.0, reference_fugacity=z0,
            energy_perturbation=1e-9)
        Random.seed!(90311)
        dO, lsO, pO = grand_canonical_nested_sampling(lsO, gc, Int64(60),
            moves, rf_save("o"))
        rf_cleanup("o")
        lsE = rf_fresh()
        ig = IdealGasReferencedGCNSParameters(mc_steps=30,
            reference_fugacity=z0, energy_perturbation=1e-9)
        Random.seed!(90311)
        dE, lsE, pE = ideal_gas_referenced_nested_sampling(lsE, ig, Int64(60),
            moves, rf_save("e"))
        rf_cleanup("e")
        @test nrow(dO) > 0
        @test nrow(dO) == nrow(dE)
        @test dO.iter == dE.iter
        @test dO.energy == dE.emax
        @test dO.num_particles == dE.num_particles
        @test dO.omega == dO.energy
        @test rf_live_E(lsO) == rf_live_E(lsE)
        @test rf_live_N(lsO) == rf_live_N(lsE)
        for k in rf_keys
            @test pO.move_stats[k] == pE.move_stats[k]
        end
    end

    @testset "exact enumeration at z0 = 0.5" begin
        # Fixture and gate of the calibrated precedent ("Validation against
        # exact grand-canonical enumeration" in test-grand_canonical_ns.jl):
        # the 4x4 two-shell lattice gas at mu = -0.05 eV and T = 300 K,
        # K = 100 walkers, 3000 iterations, mc_steps = 100, gated at rtol 0.3
        # on <E> and <N>. Here the ladder is sampled under the prior at
        # z0 = 0.5 (one site in three occupied initially, z0-weighted
        # insert/delete acceptances) and reduced with the matching keyword,
        # so the exact averages at the run's own mu must come back. The gate
        # is the precedent's own, calibrated on the uniform prior; it is on
        # the loose side by construction, and this testset had not been
        # executed when it was written. It discriminates: dropping the
        # z0^(-N) factor evaluates the ensemble at mu + kT ln z0 = -0.068 eV,
        # where the exact averages (about <N> = 7.2 and <E> = -0.38 eV against
        # 11.7 and -0.69 eV at mu) sit 38 % and 45 % away, outside the gate.
        # energy_perturbation is 1e-9 rather than the precedent's 1e-12, the
        # documented choice above the K^2 eps(E_bound) tie-breaking floor
        # (beta times 1e-9 is 4e-8, far inside the gate).
        M = 16
        ham2 = rf_ham2()
        mu = -0.05
        kb = 8.617333262e-5
        T = 300.0
        beta = 1.0 / (kb * T)

        # Exact reference: a plain loop over all 2^16 occupancies of one
        # probe lattice (no library enumeration, no RNG)
        exact_z = 0.0
        exact_E = 0.0
        exact_N = 0.0
        probe = rf_lat2()
        for mask in 0:(2^M - 1)
            for site in 1:M
                probe.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            E_val = interacting_energy(probe, ham2).val
            N_val = sum(probe.components[1])
            boltz = exp(-beta * (E_val - mu * N_val))
            exact_z += boltz
            exact_E += boltz * E_val
            exact_N += boltz * N_val
        end
        exact_mean_E = exact_E / exact_z
        exact_mean_N = exact_N / exact_z

        z0 = 0.5
        K = 100
        ls = LatticeGasWalkers(replicate_walkers(rf_lat2(), K), ham2;
                               assign_energy=false)
        gc = GrandCanonicalNestedSamplingParameters(mc_steps=100,
            chemical_potential=mu, reference_fugacity=z0,
            energy_perturbation=1e-9)
        Random.seed!(90321)
        df, _, _ = grand_canonical_nested_sampling(ls, gc, Int64(3000),
            MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25), rf_save("x"))
        rf_cleanup("x")
        @test nrow(df) > 0
        @test names(df) == ["iter", "omega", "energy", "num_particles"]
        @test issorted(df.omega, rev=true)

        mean_E_ns, Cv_ns, mean_N_ns = gc_thermodynamic_stats(
            df, [beta], K, mu; reference_fugacity=z0)
        @test isapprox(mean_E_ns[1], exact_mean_E; rtol=0.3)
        @test isapprox(mean_N_ns[1], exact_mean_N; rtol=0.3)
        @test isfinite(Cv_ns[1])
        # The reduction without the factor sits below in <N>, exactly rather
        # than statistically: the factor 2^N_j reweights toward the high-N
        # shells, and the covariance of N with an increasing function of N
        # under positive weights is positive whenever N varies over the
        # ledger, which a ladder from the one-third-occupied prior to the
        # dense ground state guarantees
        _, _, mean_N_u = gc_thermodynamic_stats(df, [beta], K, mu)
        @test length(unique(df.num_particles)) > 1
        @test mean_N_u[1] < mean_N_ns[1]
    end

    @testset "z0 = 1 reduction identity and the shell factor" begin
        # On a ledger of the Omega-sorted driver: reference_fugacity = 1.0
        # against the call that never mentions it, field for field, on both
        # methods and both compression conventions
        K = 16
        mu = -0.05
        kb = 8.617333262e-5
        ls = LatticeGasWalkers(replicate_walkers(rf_lat2(), K), rf_ham2();
                               assign_energy=false)
        gc = GrandCanonicalNestedSamplingParameters(mc_steps=30,
            chemical_potential=mu, energy_perturbation=1e-9)
        Random.seed!(90331)
        df, _, _ = grand_canonical_nested_sampling(ls, gc, Int64(200),
            MCGrandCanonicalMoves(), rf_save("z"))
        rf_cleanup("z")
        @test nrow(df) > 0
        betas = [1.0 / (kb * 300.0), 1.0 / (kb * 600.0)]
        r0 = gc_thermodynamic_stats(df, betas, K, mu)
        r1 = gc_thermodynamic_stats(df, betas, K, mu; reference_fugacity=1.0)
        @test r0[1] == r1[1]
        @test r0[2] == r1[2]
        @test r0[3] == r1[3]
        rm0 = gc_thermodynamic_stats(df, betas, K, mu; compression=:mean)
        rm1 = gc_thermodynamic_stats(df, betas, K, mu; compression=:mean,
                                     reference_fugacity=1.0)
        @test rm0[1] == rm1[1]
        @test rm0[2] == rm1[2]
        @test rm0[3] == rm1[3]
        w = ωᵢ(df.iter, K)
        v0 = gc_thermodynamic_stats(betas[1], w, df.omega, df.energy,
                                    df.num_particles, mu)
        v1 = gc_thermodynamic_stats(betas[1], w, df.omega, df.energy,
                                    df.num_particles, mu; reference_fugacity=1.0)
        @test v0 == v1
        # the DataFrame method is the vector method per beta
        @test (r0[1][1], r0[2][1], r0[3][1]) == v0
        # a non-unit z0 forwards exactly and changes the answer
        rz = gc_thermodynamic_stats(df, betas, K, mu; reference_fugacity=0.5)
        vz = gc_thermodynamic_stats(betas[1], w, df.omega, df.energy,
                                    df.num_particles, mu; reference_fugacity=0.5)
        @test (rz[1][1], rz[2][1], rz[3][1]) == vz
        @test rz[3] != r0[3]
        # positivity on both methods
        @test_throws ArgumentError gc_thermodynamic_stats(df, betas, K, mu;
                                                          reference_fugacity=0.0)
        @test_throws ArgumentError gc_thermodynamic_stats(betas[1], w, df.omega,
            df.energy, df.num_particles, mu; reference_fugacity=-1.0)
        # A two-shell ledger with a closed form: weights 0.5 each, shells
        # (E, N) = (0, 0) and (-1, 1), mu = -0.5, beta = 1, kb = 1. At
        # z0 = 0.5 the second shell's weight gains z0^(-1) = 2, so its
        # normalized weight is p = 2 e^{1/2} / (1 + 2 e^{1/2}); <N> = p,
        # <E> = -p, <E^2> = p, <EN> = -p, and
        # C = Var(E) - mu Cov(E, N) = (p - p^2) - mu (-p + p^2)
        μh = -0.5
        ωh = [0.5, 0.5]
        Ωh = [0.0, -1.0 - μh * 1.0]
        Eh = [0.0, -1.0]
        Nh = [0, 1]
        u_h, cv_h, n_h = gc_thermodynamic_stats(1.0, ωh, Ωh, Eh, Nh, μh;
                                                kb=1.0, reference_fugacity=0.5)
        p2 = 2.0 * exp(0.5) / (1.0 + 2.0 * exp(0.5))
        @test isapprox(n_h, p2; rtol=1e-12)
        @test isapprox(u_h, -p2; rtol=1e-12)
        @test isapprox(cv_h, (p2 - p2^2) - μh * (-p2 + p2^2); rtol=1e-12)
        # the same ledger under the uniform prior, and the default identity
        u_1, cv_1, n_1 = gc_thermodynamic_stats(1.0, ωh, Ωh, Eh, Nh, μh; kb=1.0)
        p1 = exp(0.5) / (1.0 + exp(0.5))
        @test isapprox(n_1, p1; rtol=1e-12)
        @test (u_1, cv_1, n_1) == gc_thermodynamic_stats(1.0, ωh, Ωh, Eh, Nh,
            μh; kb=1.0, reference_fugacity=1.0)
    end
end
