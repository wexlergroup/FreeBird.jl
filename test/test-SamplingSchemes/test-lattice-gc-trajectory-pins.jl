@testset "Lattice GC trajectory pins" begin
    using Random
    using Unitful
    using DataFrames

    # Absolute seeded bit-identity pins on the lattice grand-canonical
    # trajectories, recorded on the shipped kernels and drivers. Every other
    # seeded lattice trajectory test is relative (same-seed self- or A/B
    # equality), so none of them can falsify a "trajectories are bit-for-bit
    # unchanged" promise from kernel-internals work; these pins can.
    #
    # Determinism across CI architectures: every recorded float is produced by
    # a short, source-ordered sequence of scalar operations (the fixed nested
    # accumulation of lattice_interaction_energy, one perturbation product per
    # step, one mu*N subtraction for Omega), with no vectorized reductions, so
    # exact == pins are expected to hold on every CI leg. Disclosed fallback:
    # should any leg falsify that argument on the energy digits, those specific
    # pins drop to the atomistic rtol 1e-12 convention while every integer pin
    # stays exact.
    #
    # move_stats pins are asserted per key, by name, never as an exact key set
    # or exact-tuple comparison, so later changes that append subset counters
    # extend the key set without touching this file.
    #
    # Captured at dev 60cec59e; capture reproduced identically across two
    # separate Julia processes before recording.

    pin_lattice() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false), cutoff_radii=[1.1],
        components=[[false for _ in 1:16]], adsorptions=:full)

    pin_ham() = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")

    pin_save(tag) = SaveEveryN("t_pin_$(tag).csv", "t_pin_$(tag).traj",
                               "t_pin_$(tag).ls", 1000000, 1000000, 1000000)
    pin_cleanup(tag) = rm.(["t_pin_$(tag).csv", "t_pin_$(tag).traj",
                            "t_pin_$(tag).ls"], force=true)

    pin_keys = (:swap_attempted, :swap_accepted,
        :cluster_attempted, :cluster_accepted,
        :insert_uniform_attempted, :insert_uniform_accepted,
        :insert_biased_attempted, :insert_biased_accepted,
        :delete_attempted, :delete_accepted)

    function pin_run_igref(seed, routine, tag)
        Random.seed!(seed)
        lat = pin_lattice()
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                   for _ in 1:10]
        ls = LatticeGasWalkers(walkers, pin_ham(); assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=30,
            reference_fugacity=1.0, energy_perturbation=1e-9)
        df, _, pout = ideal_gas_referenced_nested_sampling(ls, params,
            Int64(50), routine, pin_save(tag))
        pin_cleanup(tag)
        return df, pout
    end

    function pin_run_omega(seed, tag)
        Random.seed!(seed)
        lat = pin_lattice()
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                   for _ in 1:10]
        ls = LatticeGasWalkers(walkers, pin_ham(); assign_energy=false)
        gc = GrandCanonicalNestedSamplingParameters(mc_steps=30,
            chemical_potential=-0.05, energy_perturbation=1e-9)
        df, _, pout = grand_canonical_nested_sampling(ls, gc, Int64(50),
            MCGrandCanonicalMoves(), pin_save(tag))
        pin_cleanup(tag)
        return df, pout
    end

    @testset "ideal-gas-referenced driver pin" begin
        df, pout = pin_run_igref(90231, MCGrandCanonicalMoves(), "a")
        @test nrow(df) == 50
        @test df.iter == collect(1:50)
        @test df.emax == [-0.2699999999644116, -0.35000000042328644, -0.36999999979314385, -0.39000000005592456, -0.39000000046720606, -0.41000000027337313, -0.410000000342582, -0.4199999999512617, -0.4400000004303773, -0.449999999637877, -0.4500000000943816, -0.46000000016108417, -0.4600000003650343, -0.460000000497379, -0.46999999970808093, -0.5100000000263559, -0.519999999625387, -0.5199999998185273, -0.5299999996474583, -0.5300000001131523, -0.5300000002052305, -0.5599999998829622, -0.5799999996419265, -0.580000000314528, -0.5800000003592434, -0.5899999996467072, -0.5900000004915871, -0.5999999998134353, -0.6000000002701145, -0.6499999995846582, -0.6599999995732023, -0.6599999996888175, -0.6600000002741281, -0.6600000004280429, -0.6699999995186419, -0.6699999997320831, -0.6699999998576416, -0.6699999999810813, -0.67000000014307, -0.7199999995272947, -0.719999999583073, -0.7199999996755059, -0.7200000002702595, -0.7200000003391624, -0.7200000004629845, -0.7200000004754998, -0.729999999707986, -0.7299999998527085, -0.7299999998696584, -0.7299999998737972]
        @test df.num_particles == [6, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
        expected = Dict(:swap_attempted => 741, :swap_accepted => 702,
            :cluster_attempted => 0, :cluster_accepted => 0,
            :insert_uniform_attempted => 378, :insert_uniform_accepted => 124,
            :insert_biased_attempted => 0, :insert_biased_accepted => 0,
            :delete_attempted => 381, :delete_accepted => 136)
        for k in pin_keys
            @test pout.move_stats[k] == expected[k]
        end
    end

    @testset "ideal-gas-referenced driver pin, cluster branch" begin
        df, pout = pin_run_igref(90232,
            MCGrandCanonicalMoves(clusters_freq=2, swaps_freq=2), "b")
        @test nrow(df) == 50
        @test df.iter == collect(1:50)
        @test df.emax == [-0.2300000004182935, -0.3399999997672429, -0.38999999987374245, -0.3900000004700098, -0.39999999957200477, -0.42000000018662215, -0.4200000002205665, -0.4499999995620119, -0.44999999972912813, -0.4500000001240203, -0.4600000002612164, -0.5000000002334599, -0.5100000002560634, -0.5199999999750904, -0.5200000001304839, -0.5299999999059581, -0.5300000001117027, -0.5300000002035519, -0.5300000002146336, -0.539999999686027, -0.5700000003883858, -0.5900000002727569, -0.5900000003044019, -0.5900000004969573, -0.5999999999497403, -0.6099999998167888, -0.6499999997268328, -0.6500000004233193, -0.6599999995451575, -0.6599999995611678, -0.6599999998056706, -0.660000000043681, -0.6600000001236203, -0.660000000275807, -0.6600000003731324, -0.6600000004323584, -0.660000000471613, -0.669999999508748, -0.6699999996064473, -0.6699999996878839, -0.7199999998769969, -0.7299999996237462, -0.7299999997991685, -0.7299999999686371, -0.7300000003368419, -0.7300000003425152, -0.7399999995748048, -0.7400000000570205, -0.7400000001710221, -0.7400000002212038]
        @test df.num_particles == [5, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
        expected = Dict(:swap_attempted => 405, :swap_accepted => 360,
            :cluster_attempted => 347, :cluster_accepted => 308,
            :insert_uniform_attempted => 375, :insert_uniform_accepted => 126,
            :insert_biased_attempted => 0, :insert_biased_accepted => 0,
            :delete_attempted => 373, :delete_accepted => 121)
        for k in pin_keys
            @test pout.move_stats[k] == expected[k]
        end
    end

    @testset "Omega-sorted driver pin" begin
        df, pout = pin_run_omega(90233, "c")
        @test nrow(df) == 50
        @test df.iter == collect(1:50)
        @test df.omega == [0.02000000013447678, 0.01999999998702584, 0.019999999986485273, 0.019999999984942507, 0.019999999843285043, 0.010000000234017203, 0.009999999923124003, 4.838131562046044e-10, 2.9381996835553537e-11, 2.1772195157865326e-11, -2.8480401370600816e-10, -0.009999999627352818, -0.009999999630286083, -0.009999999772195456, -0.009999999962406192, -0.010000000299792533, -0.01999999960906279, -0.01999999996136026, -0.02000000043075012, -0.029999999790141396, -0.030000000428883733, -0.04000000032482687, -0.04000000038598894, -0.04000000042551077, -0.040000000488346954, -0.04999999994894666, -0.05000000000658367, -0.05000000005244898, -0.05000000033215923, -0.059999999508962176, -0.05999999969670866, -0.059999999806152005, -0.06000000022808638, -0.06000000026363628, -0.06000000036172204, -0.06000000041743647, -0.0600000004348068, -0.06000000046202303, -0.06000000049374843, -0.0699999997062476, -0.06999999980868166, -0.06999999986918726, -0.06999999991370986, -0.06999999999901241, -0.0700000000474903, -0.07000000019467689, -0.07000000025162045, -0.07000000025832032, -0.07000000030623177, -0.07999999995250584]
        @test df.energy == [-0.27999999986552326, -0.3300000000129742, -0.33000000001351476, -0.1800000000150575, -0.280000000156715, -0.4399999997659828, -0.34000000007687603, -0.44999999951618685, -0.399999999970618, -0.34999999997822784, -0.40000000028480404, -0.45999999962735283, -0.5099999996302861, -0.3599999997721955, -0.5599999999624062, -0.41000000029979256, -0.3699999996090628, -0.5199999999613603, -0.42000000043075014, -0.4799999997901414, -0.5800000004288838, -0.5400000003248269, -0.590000000385989, -0.5400000004255108, -0.540000000488347, -0.5999999999489467, -0.6500000000065838, -0.600000000052449, -0.6500000003321593, -0.6599999995089623, -0.6599999996967088, -0.609999999806152, -0.6600000002280865, -0.6100000002636363, -0.6600000003617221, -0.6100000004174365, -0.6600000004348069, -0.6600000004620231, -0.6600000004937485, -0.6699999997062477, -0.6699999998086817, -0.6699999998691873, -0.66999999991371, -0.6699999999990125, -0.7200000000474903, -0.7200000001946769, -0.6700000002516205, -0.6700000002583204, -0.6700000003062319, -0.7299999999525059]
        @test df.num_particles == [6, 7, 7, 4, 6, 9, 7, 9, 8, 7, 8, 9, 10, 7, 11, 8, 7, 10, 8, 9, 11, 10, 11, 10, 10, 11, 12, 11, 12, 12, 12, 11, 12, 11, 12, 11, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 12, 12, 12, 13]
        expected = Dict(:swap_attempted => 749, :swap_accepted => 600,
            :cluster_attempted => 0, :cluster_accepted => 0,
            :insert_uniform_attempted => 379, :insert_uniform_accepted => 147,
            :insert_biased_attempted => 0, :insert_biased_accepted => 0,
            :delete_attempted => 372, :delete_accepted => 147)
        for k in pin_keys
            @test pout.move_stats[k] == expected[k]
        end
    end

    @testset "canonical kernel pins" begin
        # The canonical-lattice walk kernels are default-path rewrites too;
        # these pins (captured on the shipped code, reproduced identically
        # across two Julia processes) make a stream regression in either
        # kernel visible to this gate. The random-walk ceiling is binding
        # (rate 0.945: eleven rejects exercise the revert path inside the
        # pinned stream).
        Random.seed!(90235)
        lat_rw = pin_lattice()
        for i in 1:16
            lat_rw.components[1][i] = rand() < 0.5
        end
        wk_rw = LatticeWalker(lat_rw,
            energy=interacting_energy(lat_rw, pin_ham()), iter=0)
        a_rw, r_rw, _ = MC_random_walk!(200, wk_rw, pin_ham(), -0.44;
                                        energy_perturb=1e-9)
        @test a_rw == true
        @test r_rw == 0.945
        @test wk_rw.energy.val == -0.48000000037701485
        @test sum(wk_rw.configuration.components[1]) == 9

        Random.seed!(90236)
        lat_cl = pin_lattice()
        for i in 1:16
            lat_cl.components[1][i] = rand() < 0.5
        end
        wk_cl = LatticeWalker(lat_cl,
            energy=interacting_energy(lat_cl, pin_ham()), iter=0)
        a_cl, r_cl, _ = MC_cluster_walk!(100, wk_cl, pin_ham(), -0.56, 0.3;
                                         energy_perturb=1e-9)
        @test a_cl == true
        @test r_cl == 1.0
        @test wk_cl.energy.val == -0.580000000391395
        @test sum(wk_cl.configuration.components[1]) == 11
    end

    @testset "kernel-level pin" begin
        # A permissive ceiling with every channel active: swaps, clusters,
        # uniform and biased insertion, deletion. The returned rate is a
        # single division of two integers, one rounding, deterministic.
        Random.seed!(90234)
        lat = pin_lattice()
        for i in 1:16
            lat.components[1][i] = rand() < 0.5
        end
        wk = LatticeWalker(lat, energy=interacting_energy(lat, pin_ham()),
                           iter=0)
        accept, rate, wk2, cl_acc, cl_tot, ms = MC_grand_canonical_walk!(
            200, wk, pin_ham(), 1000.0, 0.0;
            p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
            clusters_freq=2, swaps_freq=2, cluster_p=0.3,
            p_bias=0.4, bias_predicate=:contact, bias_shells=1)
        @test accept == true
        @test rate == 0.88
        @test wk2.energy.val == -0.39000000036549615
        @test sum(wk2.configuration.components[1]) == 8
        @test cl_acc == 41
        @test cl_tot == 41
        expected = Dict(:swap_attempted => 41, :swap_accepted => 41,
            :cluster_attempted => 41, :cluster_accepted => 41,
            :insert_uniform_attempted => 33, :insert_uniform_accepted => 28,
            :insert_biased_attempted => 22, :insert_biased_accepted => 20,
            :delete_attempted => 63, :delete_accepted => 46)
        for k in pin_keys
            @test ms[k] == expected[k]
        end
    end
end
