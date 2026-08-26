# Absolute seeded pins for the paths the driver pins of
# test-atomistic-igref-driver-pins.jl leave uncovered, captured at dev 8bdd7556
# (the shipped SHA at the close of the bounded-support and surface rounds) and
# reproduced digit-identically across two separate Julia processes before
# recording: the SURFACE-aware ideal-gas-referenced step loop, the plain step
# loop with the Galilean burst enabled, and a DRIVER-level run entered through
# initialize=false with an observable and a dead-point callback, so the ledger
# assembly, the observable pre-sort, and the pairing check are pinned without
# routing through the reference-law initializer. Those paths had same-process
# A/B coverage only (two runs of one seed compared inside one process), so a
# stream-neutral change touching their comparator or ceiling could pass every
# shipped testset while shifting the stream.
#
# Julia-version policy (observed on the PR's CI legs at filing, which run
# Julia {lts, 1, pre} x {ubuntu-x64, macos-aarch64}): on Julia 1.10 every
# vector reproduces digit-identically on the x64 lts leg and on the aarch64
# development machine, so architecture is not what moves these fixtures. Julia
# 1.12 flips one accept/reject decision after 10 to 13 steps in each of the
# three fixtures (the same flipped trajectory on x64 and aarch64) and the
# pre-release leg flips fixture G, so the trajectory vectors (num_particles,
# emax, the sorted live-set energies, the Galilean counters) are asserted under
# VERSION < v"1.11" only. On every version the file asserts the content the
# flips leave unchanged: the iteration sequences, nonincreasing emax, the
# compression charges as same-process logs of integer ratios (the plateau entry
# and exit steps of fixtures S and D are the same on every leg), the ledger
# schema and welds, and a same-process replay identity of each fixture (two
# runs of one seed, every vector equal under ==). Within a version the vectors
# follow the scoping of the driver pins: live sets are built deterministically
# (fixed particle counts, positions from raw Random-stdlib uniforms; the
# Distributions Poisson initializer is never drawn), the pinned stream consumes
# only rand()/rand(1:n) draws, emax and live-set energies are smooth
# Lennard-Jones accumulations pinned elementwise at rtol 1e-12 (a real stream
# change moves them by orders of magnitude), num_particles and iteration
# sequences are exact integers, and the Galilean counters are exact integers
# accumulated from the burst's segment outcomes.
#
# Fixture S (seed 424271): the 10 x 10 x 15 A slab cell of the surface-route
# testsets, four frozen H, K = 12 with two empty walkers, 60 steps (the descent
# traverses the exact empty-configuration plateau: six eviction charges).
# Fixture G (seed 424272): the 12 A periodic box, K = 12 with counts 1..12,
# galilean_steps = 2, 80 steps, no plateau. Fixture D (seed 424273): the
# driver at K = 12 with counts 1..6 doubled, cutoff 2.0 sigma (inside the half
# cell, so the minimum-image warning stays silent), 60 steps.
@testset "atomistic igref surface, Galilean, and driver-path pins (absolute, captured at 8bdd7556)" begin
    using Random
    pin3_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
    pin3_pbc = (true, true, true)
    pin3_seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], pin3_box, pin3_pbc))
    pin3_mkempty() = FastSystem(cell_vectors(pin3_seed_at), periodicity(pin3_seed_at),
                                empty(position(pin3_seed_at, :)), empty(species(pin3_seed_at, :)),
                                empty(mass(pin3_seed_at, :)))
    pin3_V = 1728.0
    pin3_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
    # driver fixture: cutoff 2.0 sigma = 5 A stays inside the half cell (6 A), so the
    # driver-level minimum-image warning never fires
    pin3_lj_d = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.0)

    function pin3_liveset(counts; lj=pin3_lj)
        walkers = AtomWalker{1}[]
        for n in counts
            w = AtomWalker{1}(pin3_mkempty())
            for _ in 1:n
                pos = SVector(rand() * 12.0, rand() * 12.0, rand() * 12.0)u"Å"
                FreeBird.AbstractWalkers.insert_particle!(w, pos, :Ar)
            end
            push!(walkers, w)
        end
        return GenericAtomWalkers(walkers, lj)
    end

    # surface fixture: the 10 x 10 x 15 A slab cell of the surface-route testsets
    pin3_sbox = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
    pin3_spbc = (true, true, false)
    pin3_sV = 1500.0
    pin3_surf_sys = FastSystem(atomic_system(
        [:H => [2.5, 2.5, 2.0]u"Å", :H => [7.5, 2.5, 2.0]u"Å",
         :H => [2.5, 7.5, 2.0]u"Å", :H => [7.5, 7.5, 2.0]u"Å"], pin3_sbox, pin3_spbc))
    pin3_mksurf() = AtomWalker(deepcopy(pin3_surf_sys); freeze_species=[:H])
    pin3_smkempty() = FastSystem(cell_vectors(pin3_surf_sys), periodicity(pin3_surf_sys),
                                 empty(position(pin3_surf_sys, :)), empty(species(pin3_surf_sys, :)),
                                 empty(mass(pin3_surf_sys, :)))
    pin3_cps = CompositeParameterSets(2, [LJParameters(epsilon=0.001, sigma=2.5, cutoff=1.8, shift=true),
                                          LJParameters(epsilon=0.003, sigma=2.5, cutoff=1.8, shift=true),
                                          LJParameters(epsilon=0.01, sigma=2.5, cutoff=1.8, shift=true)])

    function pin3_surface_liveset(counts)
        walkers = AtomWalker{1}[]
        for n in counts
            w = AtomWalker{1}(pin3_smkempty())
            for _ in 1:n
                pos = SVector(rand() * 10.0, rand() * 10.0, rand() * 15.0)u"Å"
                FreeBird.AbstractWalkers.insert_particle!(w, pos, :H)
            end
            push!(walkers, w)
        end
        return LJSurfaceWalkers(walkers, pin3_cps, pin3_mksurf(); assign_energy=true)
    end

    function pin3_step_run(seed::Int, ls, z0V::Float64, mc_steps::Int, n_steps::Int,
                           species::Symbol, routine)
        Random.seed!(seed)
        params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
            reference_activity=(z0V / (species == :H ? pin3_sV : pin3_V))u"Å^-3", species=species,
            allowed_fail_count=100_000)
        iters = Int[]
        emaxs = Float64[]
        npars = Int[]
        logts = Float64[]
        for k in 1:n_steps
            iter, emax, n_par, ls, params, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(
                ls, params, routine; ns_iteration=k, z0V=z0V)
            if !(iter isa Missing)
                push!(iters, iter)
                push!(emaxs, ustrip(u"eV", emax))
                push!(npars, n_par)
                push!(logts, log_t)
            end
        end
        live = sort([ustrip(u"eV", w.energy) for w in ls.walkers])
        return iters, emaxs, npars, logts, live, params
    end

    function pin3_driver_run(seed::Int, counts, z0V::Float64, mc_steps::Int, n_steps::Int)
        Random.seed!(seed)
        ls = pin3_liveset(counts; lj=pin3_lj_d)
        params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
            reference_activity=(z0V / pin3_V)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000)
        save = SaveEveryN(df_filename="_igpin3_d.csv", wk_filename="_igpin3_d.traj.extxyz",
                          ls_filename="_igpin3_d.ls.extxyz", n_traj=10^7, n_snap=10^7, n_info=10^7)
        seen = Int[]
        df, lso, pout = ideal_gas_referenced_nested_sampling(
            ls, params, n_steps, MCAtomGrandCanonicalMoves(), save;
            observables=[:n_obs => cfg -> Float64(length(cfg))],
            dead_point_callback=(iter, walker) -> push!(seen, walker.list_num_par[1]),
            initialize=false)
        for f in ["_igpin3_d.csv", "_igpin3_d.traj.extxyz", "_igpin3_d.ls.extxyz"]
            rm(f, force=true)
        end
        live = sort([ustrip(u"eV", w.energy) for w in lso.walkers])
        return df, seen, live
    end

    pin3_counts_s = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6]
    pin3_counts_g = collect(1:12)
    pin3_counts_d = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

    PIN_S_LC_NUM = [12, 12, 12, 12, 12, 12, 12, 12, 11, 10, 9, 8, 7, 6, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_S_LC_DEN = [13, 13, 13, 13, 13, 13, 13, 13, 12, 11, 10, 9, 8, 7, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_S_EMAX = [
        4419.129797462226, 1010.5537403793226, 85.63218407922264, 3.0163181736208453,
        2.3476913169078686, 0.11304809321577121, 0.05878582739650565, 0.002627620932742753,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -5.441923076881155e-5, -6.871164967922444e-5,
        -0.0003406155365023688, -0.0004956350901303653, -0.0005382278315483823, -0.0005499439184449129,
        -0.0008148414203291925, -0.0023859551644023506, -0.0025833274219798797, -0.002730482668209874,
        -0.002861228693387175, -0.0038088447900542386, -0.0038353622008285146, -0.003910989504662127,
        -0.003928601246190308, -0.0039793192138034405, -0.004065028055560342, -0.004068628577770526,
        -0.004132363973788921, -0.004387253004923683, -0.004404960263507768, -0.004519053612516182,
        -0.0045475355557148545, -0.004618143320768944, -0.004679206712934623, -0.005275838551975641,
        -0.005523001474009049, -0.005661567362651757, -0.005835726044533011, -0.005886186642191259,
        -0.006316594270489362, -0.00659960352876567, -0.007305007236287267, -0.00733737908859124,
        -0.007356128212885969, -0.007392573418900009, -0.007412815636501858, -0.007456064595922726,
        -0.007551570801007557, -0.007577477907205688, -0.00760547936811748, -0.007869436107156166,
        -0.007997165956837693, -0.008368042409748571, -0.008386160148814305, -0.008605823929848922,
    ]
    PIN_S_NPAR = [5, 3, 6, 3, 2, 5, 5, 2, 0, 0, 1, 2, 0, 0, 4, 2, 2, 4, 1, 3, 4, 1, 2, 2, 6, 4, 3, 4, 5, 4, 4, 6, 2, 3, 4, 4, 2, 3, 3, 5, 6, 5, 4, 4, 8, 4, 7, 5, 3, 5, 4, 5, 7, 3, 4, 4, 5, 4, 5, 4]
    PIN_S_LIVE = [
        -0.01273728915228397, -0.012397731251003045, -0.010812381037838537, -0.01042208585836983,
        -0.010201130479888194, -0.010187948491559028, -0.010042644350024199, -0.009845398115053144,
        -0.00903952910411551, -0.008919187994820547, -0.008863258219377767, -0.008855509136200321,
    ]
    PIN_G_LC_NUM = [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_G_LC_DEN = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_G_EMAX = [
        3144.2884568710397, 2302.4143293077777, 392.5971503756005, 115.29551252503245,
        1.1140972574803694, 0.21647293549669352, 0.20634030197749, 0.12544620921434657,
        0.03191805811958175, 0.0, -1.0405478849798388e-5, -0.000711188684746097,
        -0.0035533214781051157, -0.004795511737167576, -0.004814838855556114, -0.005242624365113564,
        -0.005787482246932602, -0.006133187535164323, -0.008596611938561198, -0.01258731466310189,
        -0.013218932602700219, -0.015513513622867917, -0.020059509308067432, -0.02426690897259415,
        -0.02475790706388038, -0.031037054452576568, -0.03217565968972192, -0.034192245985379295,
        -0.037648036595296266, -0.03786350996042457, -0.04249446372871251, -0.04383553852707076,
        -0.04740597612998426, -0.04772182804126908, -0.048592096206524264, -0.04996145323480243,
        -0.05248632251053338, -0.053911225852430494, -0.05725050234709659, -0.05774939308789473,
        -0.0577587381231814, -0.06268916944763968, -0.06452808170782498, -0.06652815979244109,
        -0.06666585161212202, -0.06848562695298296, -0.06871125543920967, -0.0709576345379536,
        -0.07121060025500044, -0.07174746949392585, -0.07233952329767018, -0.07289241505201287,
        -0.07615416064121418, -0.07641064148494951, -0.07904787190406168, -0.07933428695796324,
        -0.07959637166385172, -0.08055405598336417, -0.0815415774411034, -0.08155880396102999,
        -0.08589333463247362, -0.08650752766420201, -0.08850080877962133, -0.08906807135780583,
        -0.0894783158964076, -0.09236073486731246, -0.09450805987317386, -0.09577578262897105,
        -0.09846266299325436, -0.10128499772779555, -0.10191221395572513, -0.1022010913403918,
        -0.10247972832643418, -0.10623961178250028, -0.10713418261444885, -0.11124146877709745,
        -0.11274791611074315, -0.11355742382232921, -0.11428135674495542, -0.11439020714004836,
    ]
    PIN_G_NPAR = [11, 10, 8, 12, 10, 7, 7, 7, 8, 1, 2, 3, 6, 3, 5, 5, 4, 8, 7, 5, 8, 6, 6, 9, 8, 9, 7, 9, 9, 8, 12, 9, 12, 10, 7, 7, 9, 14, 11, 10, 10, 14, 10, 12, 13, 14, 12, 13, 10, 12, 13, 13, 14, 14, 11, 11, 11, 15, 13, 14, 14, 16, 14, 17, 17, 11, 14, 15, 15, 14, 14, 13, 16, 14, 16, 15, 18, 13, 14, 16]
    PIN_G_LIVE = [
        -0.1592523598599889, -0.15022961846025504, -0.14055521156268097, -0.1334233199212144,
        -0.12931001478363915, -0.12729476378264942, -0.12392331361457364, -0.1224483880701329,
        -0.1221661065068472, -0.12027109753967548, -0.11897824010905407, -0.11645089682921181,
    ]
    PIN_G_STATS = (galilean_attempted = 640, galilean_accepted = 619, galilean_reflect_attempted = 121, galilean_reflect_evals = 121, galilean_reflect_accepted = 100)
    PIN_D_LC_NUM = [12, 12, 11, 10, 9, 8, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_D_LC_DEN = [13, 13, 12, 11, 10, 9, 8, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_D_EMAX = [
        0.5683032610215084, 0.04085402446129616, 0.0, 0.0,
        0.0, 0.0, 0.0, -0.0002415774487922488,
        -0.0004619323185189282, -0.0014413214393594253, -0.0026588419367326215, -0.0033412393551947826,
        -0.003610632599490484, -0.0038022863894361545, -0.005395378550775143, -0.006823931352110592,
        -0.006889445857491823, -0.009254923968391053, -0.009481858767731981, -0.011110679570031153,
        -0.012033701958010592, -0.012781277332849115, -0.015120124003479292, -0.015660710979635856,
        -0.015889455197580567, -0.016368455210328436, -0.016813215495566344, -0.017121019549627376,
        -0.018989232854541682, -0.02173716977503926, -0.02281788067358222, -0.023257434242204715,
        -0.023264271601147413, -0.023421062981890763, -0.024419145414096427, -0.0248688525864012,
        -0.02489509438718094, -0.02561903931868135, -0.02650766531827569, -0.0269615687328806,
        -0.028451202388155596, -0.029167500052373035, -0.029326768596927878, -0.02971251281006961,
        -0.030007476995403997, -0.030745653042293628, -0.031410690717822834, -0.0316169620840196,
        -0.03176555521577647, -0.03190513953104817, -0.0328335182977063, -0.03359010249234158,
        -0.03363846885340897, -0.0342850545199016, -0.035160225722636024, -0.0352934784406423,
        -0.036366085579350886, -0.03826528737529982, -0.03861222019746885, -0.03865195399157628,
    ]
    PIN_D_NPAR = [6, 4, 1, 1, 2, 2, 4, 3, 6, 6, 4, 3, 4, 5, 5, 6, 5, 3, 5, 7, 4, 7, 6, 6, 6, 7, 6, 9, 9, 8, 7, 6, 9, 8, 9, 8, 7, 5, 8, 8, 8, 7, 5, 6, 8, 8, 9, 9, 9, 8, 6, 6, 9, 7, 7, 8, 7, 6, 10, 8]
    PIN_D_LIVE = [
        -0.06807237923919043, -0.054164644627748615, -0.04796606078869758, -0.0446797791601345,
        -0.04457319055415465, -0.04373742360162495, -0.04257064966474083, -0.04174472546941582,
        -0.0409156328602522, -0.04087820938898697, -0.039766124292473684, -0.039181799166424246,
    ]

    # the trajectory vectors are asserted on Julia 1.10 only (header); the
    # replay identity and the charge, order, and weld asserts run on every version
    pin3_assert_trajectory = VERSION < v"1.11"

    @testset "fixture S: surface-aware step loop through the empty plateau (seed 424271)" begin
        Random.seed!(424271)
        ls_s = pin3_surface_liveset(pin3_counts_s)
        iters, emaxs, npars, logts, live, _ = pin3_step_run(424271, ls_s, 4.0, 40, 60, :H,
                                                            MCAtomGrandCanonicalMoves())
        @test iters == collect(1:60)
        @test issorted(emaxs, rev=true)
        @test logts == log.(PIN_S_LC_NUM ./ PIN_S_LC_DEN)
        if pin3_assert_trajectory
            @test npars == PIN_S_NPAR
            @test all(isapprox.(emaxs, PIN_S_EMAX; rtol=1e-12, atol=0.0))
            @test all(isapprox.(live, PIN_S_LIVE; rtol=1e-12, atol=0.0))
        end
        # same-process replay of the fixture
        Random.seed!(424271)
        ls_s2 = pin3_surface_liveset(pin3_counts_s)
        iters2, emaxs2, npars2, logts2, live2, _ = pin3_step_run(424271, ls_s2, 4.0, 40, 60, :H,
                                                                 MCAtomGrandCanonicalMoves())
        @test iters2 == iters
        @test npars2 == npars
        @test logts2 == logts
        @test emaxs2 == emaxs
        @test live2 == live
    end

    @testset "fixture G: plain step loop with the Galilean burst (seed 424272)" begin
        Random.seed!(424272)
        ls_g = pin3_liveset(pin3_counts_g)
        iters, emaxs, npars, logts, live, params = pin3_step_run(424272, ls_g, 12.0, 40, 80, :Ar,
            MCAtomGrandCanonicalMoves(galilean_steps=2, galilean_n_refresh=4, galilean_step_size=0.5))
        @test iters == collect(1:80)
        @test issorted(emaxs, rev=true)
        @test logts == log.(PIN_G_LC_NUM ./ PIN_G_LC_DEN)
        if pin3_assert_trajectory
            @test npars == PIN_G_NPAR
            @test all(isapprox.(emaxs, PIN_G_EMAX; rtol=1e-12, atol=0.0))
            @test all(isapprox.(live, PIN_G_LIVE; rtol=1e-12, atol=0.0))
            for (k, v) in pairs(PIN_G_STATS)
                @test params.move_stats[k] == v
            end
        end
        # same-process replay of the fixture
        Random.seed!(424272)
        ls_g2 = pin3_liveset(pin3_counts_g)
        iters2, emaxs2, npars2, logts2, live2, params2 = pin3_step_run(424272, ls_g2, 12.0, 40, 80, :Ar,
            MCAtomGrandCanonicalMoves(galilean_steps=2, galilean_n_refresh=4, galilean_step_size=0.5))
        @test iters2 == iters
        @test npars2 == npars
        @test logts2 == logts
        @test emaxs2 == emaxs
        @test live2 == live
        @test all(params2.move_stats[k] == params.move_stats[k] for k in keys(PIN_G_STATS))
    end

    @testset "fixture D: driver entered through initialize=false with an observable (seed 424273)" begin
        df, seen, live = pin3_driver_run(424273, pin3_counts_d, 6.0, 40, 60)
        @test nrow(df) == 60
        @test names(df) == ["iter", "emax", "num_particles", "log_compression", "n_obs"]
        @test df.iter == collect(1:60)
        @test issorted(df.emax, rev=true)
        @test df.log_compression == log.(PIN_D_LC_NUM ./ PIN_D_LC_DEN)
        if pin3_assert_trajectory
            @test df.num_particles == PIN_D_NPAR
            @test all(isapprox.(df.emax, PIN_D_EMAX; rtol=1e-12, atol=0.0))
            @test all(isapprox.(live, PIN_D_LIVE; rtol=1e-12, atol=0.0))
        end
        @test df.n_obs == Float64.(df.num_particles)
        @test seen == df.num_particles
        # same-process replay of the fixture
        df2, seen2, live2 = pin3_driver_run(424273, pin3_counts_d, 6.0, 40, 60)
        @test df2 == df
        @test seen2 == seen
        @test live2 == live
    end
end
