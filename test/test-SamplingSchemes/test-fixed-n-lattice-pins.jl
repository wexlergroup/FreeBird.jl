@testset "Fixed-N lattice trajectory pins" begin
    using Random
    using Unitful
    using DataFrames

    # Absolute seeded bit-identity pins on the fixed-N (canonical) lattice
    # walk kernel `MC_random_walk!` and on the two nested-sampling drivers
    # that reach it, `MCRandomWalkClone` and the lattice `MCMixedMoves` step,
    # recorded on the shipped code. They extend the square-lattice canonical
    # kernel pin of test-lattice-gc-trajectory-pins.jl to the cells the
    # fixed-N kernel is used on beyond the square lattice (a two-site
    # triangular cell and an aligned two-layer triangular cell with a
    # direction-resolved shell ladder) and to a multi-body Hamiltonian, so
    # that kernel-internals work on the fixed-N walk (an incremental energy
    # path, for instance) has a gate that can falsify a "the default path is
    # stream-neutral" promise: every seeded fixed-N test elsewhere is a
    # same-seed self- or A/B comparison, which such a change passes trivially.
    #
    # Determinism across CI architectures: every recorded float is produced
    # by a short, source-ordered sequence of scalar operations (the fixed
    # nested accumulation of the pair energy, the per-embedding scalar count
    # of `cluster_energy` times one coupling, one perturbation product per
    # step, plus in the drivers the stable sort of the live set), with no
    # vectorized reductions, so exact == pins are expected to hold on every
    # CI leg. The neighbor lists behind them are integer-valued cutoff
    # decisions with the nearest excluded distance at least 0.15 from every
    # cutoff. Disclosed fallback: should any leg falsify that argument on
    # the energy digits, those specific pins drop to the atomistic rtol
    # 1e-12 convention while every integer pin stays exact.
    #
    # Fixtures are built deterministically: fixed particle counts filled by
    # raw `rand` site draws (no sampler library), so the pinned stream is
    # the fixture fill plus the kernel's own draws. Ceilings are chosen so
    # that every kernel pin rejects part of its 200 proposals (rates 0.85,
    # 0.825, 0.86), exercising the revert path inside the pinned stream, and
    # the driver runs stop before the 4 x 4 ground state, where consecutive
    # failed iterations would otherwise dominate the ledger; each run loses
    # a few iterations (3, 1 and 6), so the ledger length carries the
    # failure count and `df.iter` is pinned against it.
    #
    # Captured at dev 60113fd1; capture reproduced identically across two
    # separate Julia processes before recording.

    pin_square() = MLattice{1,SquareLattice}(lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
        periodicity=(true, true, false), cutoff_radii=[1.1],
        components=[[false for _ in 1:16]], adsorptions=:full)

    # Two-site triangular cell, 48 sites, one nearest-neighbor shell of six.
    pin_tri() = MLattice{1,TriangularLattice}(lattice_constant=1.0,
        supercell_dimensions=(6, 4, 1), periodicity=(true, true, false),
        cutoff_radii=[1.1], components=[[false for _ in 1:48]],
        adsorptions=:full)

    # Aligned two-layer triangular cell, 64 sites, spacing 1.25 with the
    # ladder [1.1, 1.3]: shell 1 is the six in-plane neighbors, shell 2 the
    # single axial neighbor (the next distances, 1.6 and sqrt(3), are
    # excluded), so the two couplings are direction-resolved.
    pin_tri_layered() = MLattice{1,TriangularLattice}(lattice_constant=1.0,
        interlayer_spacing=1.25, supercell_dimensions=(4, 4, 2),
        periodicity=(true, true, false), cutoff_radii=[1.1, 1.3],
        components=[[false for _ in 1:64]], adsorptions=:full)

    pin_ham() = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
    pin_ham_layered() = GenericLatticeHamiltonian(-0.04, [-0.01, -0.015], u"eV")

    # Fixed-count fill by rejection on raw site draws; the draw count is
    # part of the pinned stream.
    function pin_fill!(lat, n)
        occ = lat.components[1]
        fill!(occ, false)
        k = 0
        while k < n
            s = rand(eachindex(occ))
            if !occ[s]
                occ[s] = true
                k += 1
            end
        end
        return lat
    end

    pin_save(tag) = SaveEveryN("t_fnpin_$(tag).csv", "t_fnpin_$(tag).traj",
                               "t_fnpin_$(tag).ls", 1000000, 1000000, 1000000)
    pin_cleanup(tag) = rm.(["t_fnpin_$(tag).csv", "t_fnpin_$(tag).traj",
                            "t_fnpin_$(tag).ls"], force=true)

    # Twenty walkers of eight particles on the 4 x 4 cell, 30 steps per
    # iteration, 120 iterations; energies assigned unperturbed (the walk
    # perturbs), so the initial live set carries no random draw beyond the
    # fill.
    function pin_run(seed, routine, tag)
        Random.seed!(seed)
        h = pin_ham()
        walkers = [LatticeWalker(pin_fill!(pin_square(), 8), energy=0.0u"eV",
                                 iter=0) for _ in 1:20]
        for w in walkers
            w.energy = interacting_energy(w.configuration, h)
        end
        ls = LatticeGasWalkers(walkers, h; assign_energy=false)
        params = LatticeNestedSamplingParameters(mc_steps=30,
                                                 energy_perturbation=1e-9)
        df, lsf, _ = nested_sampling(ls, params, Int64(120), routine,
                                     pin_save(tag))
        pin_cleanup(tag)
        return df, lsf
    end

    @testset "fixed-N driver pin, MCRandomWalkClone" begin
        df, lsf = pin_run(90251, MCRandomWalkClone(), "a")
        @test nrow(df) == 117
        @test df.iter == collect(1:117)
        @test df.emax == [-0.38, -0.38, -0.38, -0.38, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39, -0.3900000001807564, -0.39000000023000353, -0.3900000004261791, -0.3900000004711595, -0.3999999995047626, -0.39999999955447935, -0.39999999974899575, -0.39999999981163703, -0.39999999988926627, -0.3999999999137062, -0.4, -0.4, -0.4, -0.4, -0.4, -0.40000000000790964, -0.4000000000606414, -0.4000000001458669, -0.4000000001664976, -0.40000000018040094, -0.40000000018047827, -0.4000000002054431, -0.4000000002605866, -0.4000000002679616, -0.40000000033796734, -0.4000000004471777, -0.4000000004915771, -0.4099999995431646, -0.40999999955797045, -0.4099999996701965, -0.40999999981739, -0.40999999983727325, -0.4099999998660656, -0.4099999998778149, -0.41000000000000003, -0.41000000000000003, -0.4100000001544447, -0.41000000019585775, -0.4100000002134545, -0.41000000021718475, -0.41000000021931515, -0.41000000025201205, -0.4100000003117903, -0.41000000032032025, -0.4100000003870923, -0.4199999995325914, -0.41999999961832873, -0.41999999963905355, -0.4199999996439117, -0.4199999997175096, -0.4199999997763902, -0.4199999997801532, -0.4199999997990529, -0.41999999983686853, -0.4199999998780859, -0.41999999990002973, -0.41999999991617654, -0.41999999999177834, -0.420000000000515, -0.42000000000398846, -0.4200000000129841, -0.4200000000281086, -0.4200000000296632, -0.4200000000302263, -0.42000000008692856, -0.42000000009659205, -0.42000000009795174, -0.4200000001075439, -0.42000000010939464, -0.42000000011862665, -0.42000000013169486, -0.4200000001497417, -0.42000000019180794, -0.42000000023487105, -0.42000000025526546, -0.4200000002558304, -0.42000000025692175, -0.42000000027108836, -0.42000000027121315, -0.4200000002903444, -0.420000000298263, -0.42000000030389606, -0.4200000003139311, -0.4200000003271934, -0.42000000032939416, -0.4200000003326725, -0.4200000003363392, -0.42000000034024054, -0.42000000034219265, -0.4200000003486644, -0.4200000003669525, -0.4200000003744123, -0.42000000038615504, -0.42000000041661206, -0.42000000042932256, -0.4200000004332465, -0.4200000004339031, -0.4200000004393559, -0.42000000044205466, -0.4200000004471257, -0.4200000004507202, -0.4200000004512742, -0.4200000004559789, -0.42000000045981334]
        @test all(sum(w.configuration.components[1]) == 8 for w in lsf.walkers)
    end

    @testset "fixed-N driver pin, MCMixedMoves swaps only" begin
        # The lattice mixed step with clusters_freq = 0 reaches the same
        # kernel through a different step method (n_local = mc_steps).
        df, lsf = pin_run(90252, MCMixedMoves(walks_freq=1, clusters_freq=0), "b")
        @test nrow(df) == 119
        @test df.iter == collect(1:119)
        @test df.emax == [-0.37, -0.38, -0.38, -0.38, -0.38, -0.38, -0.38000000012044044, -0.39, -0.39, -0.3900000002040753, -0.3900000004155881, -0.3900000004604998, -0.399999999540736, -0.399999999541944, -0.39999999967318084, -0.39999999969988537, -0.3999999997536736, -0.3999999999983048, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4000000000143525, -0.4000000000619271, -0.40000000007123593, -0.4000000001205468, -0.4000000001441698, -0.4000000001895667, -0.40000000023343785, -0.40000000027590693, -0.4000000002815449, -0.40000000030860877, -0.4000000003387519, -0.40000000037984484, -0.40000000039924316, -0.4000000004021102, -0.4000000004181006, -0.40000000047422624, -0.4000000004745936, -0.4099999995236799, -0.40999999960032724, -0.40999999960633715, -0.40999999963647044, -0.4099999996702715, -0.4099999998653411, -0.40999999987735203, -0.40999999989962205, -0.409999999923576, -0.41000000000000003, -0.41000000000000003, -0.41000000000000003, -0.41000000000000003, -0.4100000000833317, -0.41000000012344145, -0.4100000001311785, -0.4100000001604217, -0.4100000002191982, -0.41000000023275723, -0.4100000002734969, -0.41000000032392453, -0.41000000032454076, -0.4100000003925644, -0.4100000004626971, -0.4100000004830845, -0.4100000004854553, -0.41000000048985363, -0.41999999955717743, -0.4199999995941008, -0.4199999995995564, -0.4199999996011411, -0.41999999969848273, -0.4199999997208097, -0.41999999978197244, -0.41999999985195574, -0.4199999998819259, -0.4199999998945075, -0.41999999991593395, -0.41999999991724324, -0.4199999999397287, -0.41999999994235854, -0.41999999998171533, -0.41999999999292653, -0.42000000000000004, -0.42000000004301, -0.42000000010454486, -0.42000000011727856, -0.4200000001425728, -0.42000000016359007, -0.4200000001690306, -0.42000000020006245, -0.42000000026649786, -0.42000000027585555, -0.42000000028225515, -0.4200000002881306, -0.4200000002906304, -0.42000000029951085, -0.42000000029954154, -0.42000000030473283, -0.42000000032086676, -0.4200000003357172, -0.42000000036820956, -0.42000000038672164, -0.4200000004023818, -0.4200000004110968, -0.4200000004111511, -0.42000000041569524, -0.4200000004275372, -0.4200000004283438, -0.4200000004295292, -0.42000000043109753, -0.4200000004340939, -0.42000000043932195, -0.42000000043936, -0.4200000004435015, -0.4200000004498948, -0.4200000004578165]
        @test all(sum(w.configuration.components[1]) == 8 for w in lsf.walkers)
    end

    @testset "fixed-N driver pin, MCMixedMoves with cluster moves" begin
        # Half the steps are geometric cluster moves under the adaptive
        # growth probability (it adjusts at the 50- and 100-iteration
        # windows), half are hop-pair swaps.
        df, lsf = pin_run(90253, MCMixedMoves(walks_freq=1, clusters_freq=1), "c")
        @test nrow(df) == 114
        @test df.iter == collect(1:114)
        @test df.emax == [-0.37, -0.38, -0.3899999999341637, -0.3899999999697512, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39, -0.3900000000263392, -0.39000000023920794, -0.3900000003830612, -0.3999999996727031, -0.3999999997469902, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4000000000069155, -0.4000000000222916, -0.40000000004648734, -0.4000000000613675, -0.4000000000640361, -0.40000000008046777, -0.4000000002519381, -0.4000000003349598, -0.400000000367241, -0.4000000003856155, -0.40000000038801486, -0.4000000004114562, -0.40999999955871724, -0.40999999967874423, -0.4099999996829018, -0.40999999979236984, -0.409999999878783, -0.40999999987962304, -0.40999999991758895, -0.4099999999989046, -0.41000000000000003, -0.4100000001562529, -0.4100000001988061, -0.41000000022277494, -0.410000000248584, -0.4100000002752026, -0.4100000003054737, -0.41000000040759993, -0.41000000042751145, -0.41000000048285656, -0.4199999995136147, -0.41999999953687794, -0.419999999585839, -0.41999999964183193, -0.4199999998604259, -0.41999999988156567, -0.41999999996088366, -0.41999999997170945, -0.4199999999884855, -0.42000000000000004, -0.42000000000000004, -0.4200000000113927, -0.4200000000117214, -0.42000000001583876, -0.4200000000344042, -0.4200000000693463, -0.4200000001191661, -0.42000000012304717, -0.4200000001529655, -0.4200000001989329, -0.42000000020761874, -0.42000000021347317, -0.4200000002615855, -0.42000000028113066, -0.42000000028418993, -0.4200000003073207, -0.4200000003306373, -0.42000000033401647, -0.42000000034195434, -0.4200000003477944, -0.42000000035920393, -0.42000000036967444, -0.42000000037159224, -0.42000000037475566, -0.4200000003767698, -0.42000000038162943, -0.42000000038229185, -0.4200000004019533, -0.4200000004072921, -0.4200000004097885, -0.42000000040994934, -0.4200000004135909, -0.42000000041594965, -0.4200000004200269, -0.4200000004201385, -0.4200000004203728, -0.4200000004312553, -0.42000000043269325, -0.42000000043882224, -0.42000000044083424, -0.42000000044372443, -0.42000000045791325, -0.4200000004611362, -0.4200000004648871, -0.4200000004658448, -0.4200000004716987, -0.42000000047518604, -0.42000000047612174, -0.4200000004768893, -0.4200000004866568]
        @test all(sum(w.configuration.components[1]) == 8 for w in lsf.walkers)
    end

    @testset "fixed-N kernel pins, triangular cells" begin
        # The start energy pins the seeded fill (an exact decimal from the
        # closed-form count of occupied bonds); the ceiling sits one bond
        # energy above it so the walk both accepts and rejects.
        Random.seed!(90241)
        lat_tri = pin_fill!(pin_tri(), 24)
        wk_tri = LatticeWalker(lat_tri,
            energy=interacting_energy(lat_tri, pin_ham()), iter=0)
        @test wk_tri.energy.val == -1.33
        a_tri, r_tri, _ = MC_random_walk!(200, wk_tri, pin_ham(), -1.325;
                                          energy_perturb=1e-9)
        @test a_tri == true
        @test r_tri == 0.85
        @test wk_tri.energy.val == -1.3400000003697636
        @test sum(wk_tri.configuration.components[1]) == 24

        Random.seed!(90242)
        lat_lay = pin_fill!(pin_tri_layered(), 32)
        wk_lay = LatticeWalker(lat_lay,
            energy=interacting_energy(lat_lay, pin_ham_layered()), iter=0)
        @test wk_lay.energy.val == -1.895
        a_lay, r_lay, _ = MC_random_walk!(200, wk_lay, pin_ham_layered(), -1.8875;
                                          energy_perturb=1e-9)
        @test a_lay == true
        @test r_lay == 0.825
        @test wk_lay.energy.val == -1.9400000001221305
        @test sum(wk_lay.configuration.components[1]) == 32
    end

    @testset "fixed-N kernel pin, cluster Hamiltonian" begin
        # One repulsive right-isosceles trio figure on the 4 x 4 cell (four
        # embeddings per plaquette) over the pair Hamiltonian; the cell is
        # not a faithful quotient for this figure, so the enumeration warns
        # and follows the minimum-image convention, which is the convention
        # the energy kernel then evaluates.
        sq = pin_square()
        embs = @test_logs (:warn, r"faithful quotient") match_mode=:any enumerate_motif_embeddings(
            sq, motif_distances([(0, 0), (1, 0), (0, 1)]); expected_count=64)
        @test length(embs) == 64
        hc = ClusterLatticeHamiltonian(pin_ham(),
                                       [ClusterInteraction(0.02u"eV", embs)])
        Random.seed!(90243)
        pin_fill!(sq, 8)
        wk_cl = LatticeWalker(sq, energy=interacting_energy(sq, hc), iter=0)
        @test wk_cl.energy.val == -0.25
        a_cl, r_cl, _ = MC_random_walk!(200, wk_cl, hc, -0.245;
                                        energy_perturb=1e-9)
        @test a_cl == true
        @test r_cl == 0.86
        @test wk_cl.energy.val == -0.2600000000310176
        @test sum(wk_cl.configuration.components[1]) == 8
    end
end
