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
    # separate Julia processes before recording. The grand-canonical driver and
    # kernel pins were re-recorded for #287 on dev 1130c74e on the per-trial
    # Metropolis uniform stream (every pin whose trajectory passes through an
    # insertion or deletion proposal; the canonical kernel pins are unchanged),
    # again reproduced identically across two separate Julia processes.

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
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        df, pout = pin_run_igref(90231, MCGrandCanonicalMoves(), "a")
        @test nrow(df) == 50
        @test df.iter == collect(1:50)
        @test df.emax == [-0.2699999999644116, -0.28000000007155673, -0.35000000042328644, -0.39000000005592456, -0.39000000046720606, -0.410000000342582, -0.4500000000943816, -0.45000000009646907, -0.46000000007354164, -0.46000000016108417, -0.4600000002061177, -0.4600000003631707, -0.46999999986789, -0.4699999999706489, -0.470000000149826, -0.48000000019996863, -0.5000000004800368, -0.5100000002406192, -0.5199999999736228, -0.5200000004447404, -0.5200000004694979, -0.5299999995947576, -0.5299999999531786, -0.5299999999775558, -0.5300000002052305, -0.5300000003663144, -0.5699999998778553, -0.5800000003476262, -0.5899999997255442, -0.5899999998074682, -0.5899999999223221, -0.5900000003088514, -0.6099999999291209, -0.6100000003289007, -0.6499999995868154, -0.6499999996477377, -0.650000000317488, -0.6500000004633409, -0.659999999527596, -0.659999999709715, -0.6599999998240963, -0.6599999998790677, -0.6599999999234855, -0.6699999995186419, -0.6699999996486924, -0.6699999997320831, -0.6699999998982838, -0.6800000003683231, -0.7199999995124136, -0.7199999998002508]
        @test df.num_particles == [6, 6, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13]
        expected = Dict(:swap_attempted => 769, :swap_accepted => 719,
            :cluster_attempted => 0, :cluster_accepted => 0,
            :insert_uniform_attempted => 361, :insert_uniform_accepted => 138,
            :insert_biased_attempted => 0, :insert_biased_accepted => 0,
            :delete_attempted => 370, :delete_accepted => 149)
        for k in pin_keys
            @test pout.move_stats[k] == expected[k]
        end
    end

    @testset "ideal-gas-referenced driver pin, cluster branch" begin
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        df, pout = pin_run_igref(90232,
            MCGrandCanonicalMoves(clusters_freq=2, swaps_freq=2), "b")
        @test nrow(df) == 50
        @test df.iter == collect(1:50)
        @test df.emax == [-0.2300000004182935, -0.3399999997672429, -0.38999999987374245, -0.39000000041117894, -0.3900000004930578, -0.39999999957200477, -0.399999999711398, -0.39999999975845324, -0.39999999985338935, -0.42000000018662215, -0.4200000002205665, -0.44999999972912813, -0.44999999994061424, -0.4500000001240203, -0.46000000017753423, -0.4699999995734818, -0.4699999996027105, -0.48000000029233936, -0.5099999995079432, -0.5099999996772858, -0.510000000135971, -0.5199999996398718, -0.5199999996546878, -0.5199999996685933, -0.5200000000650293, -0.5200000001304839, -0.5200000004331922, -0.5200000004741844, -0.5299999998694175, -0.5299999999123578, -0.5300000004164148, -0.5499999999345818, -0.5799999995210208, -0.5999999996449518, -0.599999999657436, -0.6099999996431783, -0.6100000000988276, -0.6500000001118634, -0.6500000002582796, -0.6500000002598414, -0.6500000004233193, -0.6599999996704778, -0.6599999998482745, -0.660000000186492, -0.660000000191554, -0.6600000002264781, -0.6699999997122547, -0.6699999998305235, -0.6699999999711298, -0.6799999997124324]
        @test df.num_particles == [5, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
        expected = Dict(:swap_attempted => 389, :swap_accepted => 366,
            :cluster_attempted => 356, :cluster_accepted => 338,
            :insert_uniform_attempted => 365, :insert_uniform_accepted => 162,
            :insert_biased_attempted => 0, :insert_biased_accepted => 0,
            :delete_attempted => 390, :delete_accepted => 173)
        for k in pin_keys
            @test pout.move_stats[k] == expected[k]
        end
    end

    @testset "Omega-sorted driver pin" begin
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        df, pout = pin_run_omega(90233, "c")
        @test nrow(df) == 50
        @test df.iter == collect(1:50)
        @test df.omega == [0.02000000013447678, 0.01999999998702584, 0.019999999986724193, 0.019999999986485273, 0.019999999984942507, 0.019999999843285043, 0.019999999774085897, 0.010000000234017203, 0.010000000085310712, 2.1772195157865326e-11, -8.499412285090102e-11, -1.3729173353738133e-10, -4.563342481667121e-10, -0.009999999627352818, -0.00999999971154375, -0.009999999790440195, -0.010000000292133493, -0.010000000434231937, -0.01999999966160615, -0.01999999998612756, -0.020000000234044024, -0.020000000292823006, -0.020000000297820952, -0.020000000387666694, -0.020000000390253625, -0.029999999808422828, -0.030000000164860374, -0.030000000232339286, -0.030000000428883733, -0.039999999608890335, -0.03999999962735756, -0.0399999996590753, -0.03999999988496816, -0.04000000020377603, -0.04000000023650807, -0.04000000028126638, -0.0400000003512343, -0.04999999984662984, -0.04999999991830806, -0.05999999964975822, -0.059999999688591044, -0.05999999973306025, -0.05999999992375549, -0.05999999999763106, -0.060000000258337405, -0.06000000031757258, -0.0600000004185548, -0.06999999955976399, -0.06999999963697778, -0.07000000028100983]
        @test df.energy == [-0.27999999986552326, -0.3300000000129742, -0.28000000001327585, -0.33000000001351476, -0.1800000000150575, -0.280000000156715, -0.3800000002259141, -0.4399999997659828, -0.28999999991468933, -0.34999999997822784, -0.30000000008499417, -0.40000000013729176, -0.45000000045633426, -0.45999999962735283, -0.45999999971154376, -0.5099999997904402, -0.4100000002921335, -0.46000000043423195, -0.5199999996616061, -0.5199999999861276, -0.47000000023404404, -0.520000000292823, -0.47000000029782096, -0.5200000003876667, -0.5200000003902536, -0.47999999980842284, -0.5300000001648604, -0.4800000002323393, -0.5800000004288838, -0.5899999996088904, -0.5899999996273576, -0.5399999996590753, -0.5899999998849682, -0.5900000002037761, -0.6400000002365082, -0.5400000002812664, -0.6400000003512344, -0.5999999998466299, -0.6499999999183081, -0.6599999996497583, -0.6599999996885911, -0.6599999997330603, -0.6099999999237555, -0.6599999999976311, -0.6600000002583375, -0.6600000003175727, -0.6600000004185549, -0.6699999995597641, -0.6699999996369779, -0.6700000002810099]
        @test df.num_particles == [6, 7, 6, 7, 4, 6, 8, 9, 6, 7, 6, 8, 9, 9, 9, 10, 8, 9, 10, 10, 9, 10, 9, 10, 10, 9, 10, 9, 11, 11, 11, 10, 11, 11, 12, 10, 12, 11, 12, 12, 12, 12, 11, 12, 12, 12, 12, 12, 12, 12]
        expected = Dict(:swap_attempted => 726, :swap_accepted => 630,
            :cluster_attempted => 0, :cluster_accepted => 0,
            :insert_uniform_attempted => 393, :insert_uniform_accepted => 178,
            :insert_biased_attempted => 0, :insert_biased_accepted => 0,
            :delete_attempted => 381, :delete_accepted => 171)
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
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
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
        @test rate == 0.89
        @test wk2.energy.val == -0.23000000004642096
        @test sum(wk2.configuration.components[1]) == 5
        @test cl_acc == 46
        @test cl_tot == 46
        expected = Dict(:swap_attempted => 51, :swap_accepted => 51,
            :cluster_attempted => 46, :cluster_accepted => 46,
            :insert_uniform_attempted => 27, :insert_uniform_accepted => 24,
            :insert_biased_attempted => 17, :insert_biased_accepted => 16,
            :delete_attempted => 59, :delete_accepted => 41)
        for k in pin_keys
            @test ms[k] == expected[k]
        end
    end
end
