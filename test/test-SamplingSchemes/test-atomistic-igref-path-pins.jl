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
# Re-recorded for #287 on dev 1130c74e, where every insertion and deletion
# proposal of the grand-canonical kernels draws its Metropolis uniform with
# the proposal whatever the ceiling outcome and the ratio: the same fixtures,
# seeds, tolerances and version policy, every vector again reproduced
# digit-identically in two separate Julia 1.10.4 processes before recording
# (the version-policy observations below were made on the previous stream).
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
            allowed_fail_count=100_000, compression=:mean)  # compression keyword: fixture on the historical mean convention (compression=:mean); the geometric default is covered by test-compression-convention.jl
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
            allowed_fail_count=100_000, compression=:mean)  # compression keyword: fixture on the historical mean convention (compression=:mean); the geometric default is covered by test-compression-convention.jl
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

    PIN_S_LC_NUM = [12, 12, 12, 12, 12, 12, 11, 10, 9, 8, 7, 6, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_S_LC_DEN = [13, 13, 13, 13, 13, 13, 12, 11, 10, 9, 8, 7, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_S_EMAX = [
        4419.129797462226, 1010.5537403793226, 85.63218407922264, 3.0163181736208453,
        0.14772194349982362, 0.002627620932742753, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        -0.0002146729776614182, -0.0004672113308036267, -0.0004956350901303653, -0.0005382278315483823,
        -0.0009116867268160558, -0.0010912657388625871, -0.0012881399964247077, -0.0017084310373231416,
        -0.0021778952218127682, -0.002301079866406946, -0.0028786955590249814, -0.0035146970010662584,
        -0.0036413049565436287, -0.003942985263373215, -0.0042345450769048915, -0.004256581148693766,
        -0.004280386985191821, -0.004511278712621556, -0.0045475355557148545, -0.0045820594731876265,
        -0.005171159246455122, -0.0051907562800847665, -0.005331737773783143, -0.005466871798406079,
        -0.005552295330090429, -0.005606905103181623, -0.005613051447877659, -0.005835726044533011,
        -0.00614498824656448, -0.006280694875627327, -0.0066450676955275975, -0.006850342280800237,
        -0.007059367601543081, -0.00709213621678249, -0.007097452709974818, -0.007156842929153998,
        -0.007432697538399974, -0.007478866129994553, -0.008506178457092978, -0.008796681586377909,
        -0.009161048699775949, -0.009188508539292825, -0.00919131126822274, -0.009582161005415422,
        -0.009910769984324976, -0.009919075852945035, -0.010314079225674969, -0.010767437600434963,
    ]
    PIN_S_NPAR = [5, 3, 6, 3, 4, 2, 0, 0, 1, 1, 1, 1, 2, 3, 4, 1, 5, 4, 4, 4, 4, 5, 5, 4, 4, 6, 2, 2, 3, 5, 2, 2, 5, 7, 4, 7, 4, 4, 3, 4, 5, 5, 2, 5, 4, 2, 4, 5, 4, 4, 5, 4, 5, 7, 8, 4, 5, 6, 5, 7]
    PIN_S_LIVE = [
        -0.013645233940264232, -0.013597069196576266, -0.012764052846093202, -0.012378631840658237,
        -0.01225148341010435, -0.01206213470098592, -0.011980699326793203, -0.011755292014167239,
        -0.011392223107848109, -0.01138674553966261, -0.011002417701374846, -0.010850251140767151,
    ]
    PIN_G_LC_NUM = [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_G_LC_DEN = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_G_EMAX = [
        3144.2884568710397, 2302.4143293077777, 115.29551252503245, 0.9775004896201893,
        0.34839992478688314, 0.21647293549669352, 0.1316701639020404, 0.03191805811958175,
        0.0, -1.0405478849798388e-5, -0.0008734012898841969, -0.0035533214781051157,
        -0.004795511737167576, -0.005242624365113564, -0.005787482246932602, -0.007824337082524734,
        -0.01087567570357035, -0.012810734983911576, -0.014090564818902099, -0.01899993662900737,
        -0.023425262146547842, -0.02426690897259415, -0.024789165964468583, -0.026594252232910173,
        -0.02828668460324274, -0.03079135103685431, -0.034445137735324885, -0.036007225272950416,
        -0.03659882911001382, -0.039382044792129894, -0.04013662476985624, -0.0405737397417517,
        -0.04066744989147483, -0.0432610826990129, -0.04375432669540134, -0.044815281702039175,
        -0.046222084412227415, -0.04676372244246161, -0.04876372109807747, -0.05292836318282909,
        -0.054266307041118796, -0.0573446605292631, -0.06109243723152222, -0.06337025738619448,
        -0.06422503218368007, -0.0654002387539533, -0.06717268391455375, -0.06816931023843044,
        -0.06834435607611919, -0.06840032170257854, -0.07094325145095987, -0.07213082720658966,
        -0.07223368742834835, -0.07297144763733401, -0.07309218763552301, -0.07323420258946636,
        -0.07364616643361556, -0.07739027878855559, -0.0781799765067724, -0.08717318791515129,
        -0.08737524035597359, -0.08743930490035719, -0.09143119482753695, -0.09161750365637465,
        -0.09594686208807711, -0.09660835157196244, -0.09757285826118989, -0.09878947264537262,
        -0.10447618800421987, -0.10894735914567126, -0.1099958408354087, -0.11110678925671631,
        -0.11266787147694896, -0.11563739194978362, -0.11694760904811853, -0.11735718091244528,
        -0.11810172602482248, -0.11837378324343768, -0.120126131455959, -0.12503144584018971,
    ]
    PIN_G_NPAR = [11, 10, 12, 10, 9, 7, 5, 8, 1, 2, 2, 6, 3, 5, 4, 6, 5, 6, 7, 9, 10, 9, 10, 8, 6, 10, 9, 10, 9, 11, 8, 9, 11, 16, 8, 10, 11, 13, 9, 12, 10, 11, 11, 11, 11, 11, 11, 14, 12, 12, 14, 13, 13, 12, 16, 13, 14, 11, 11, 13, 14, 14, 17, 13, 14, 16, 15, 13, 14, 15, 15, 19, 15, 17, 15, 14, 16, 17, 15, 14]
    PIN_G_LIVE = [
        -0.17982725854133166, -0.15733906660728755, -0.14732361879333847, -0.14124674415504077,
        -0.13779231984277532, -0.1370371528105169, -0.13622309012941924, -0.13375024533846486,
        -0.1327236794363732, -0.1293888062073271, -0.126918093014621, -0.1252437223102858,
    ]
    PIN_G_STATS = (galilean_attempted = 640, galilean_accepted = 618, galilean_reflect_attempted = 118, galilean_reflect_evals = 118, galilean_reflect_accepted = 96)
    PIN_D_LC_NUM = [12, 12, 11, 10, 9, 8, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_D_LC_DEN = [13, 13, 12, 11, 10, 9, 8, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_D_EMAX = [
        0.5683032610215084, 0.04085402446129616, 0.0, 0.0,
        0.0, 0.0, 0.0, -0.0002415774487922488,
        -0.0011511104060799291, -0.0014413214393594253, -0.002063082636245717, -0.0033412393551947826,
        -0.0038022863894361545, -0.004792967437554651, -0.006889445857491823, -0.006921179251914143,
        -0.007165943448141244, -0.007688542761035316, -0.008607833774267968, -0.00941226302973757,
        -0.009555592840985064, -0.010092251864288762, -0.010407980510568772, -0.010720055306764259,
        -0.01199983784088041, -0.012077571565124817, -0.013333345311671805, -0.01343278164860811,
        -0.013997499403684624, -0.014127108567379297, -0.014138375742291783, -0.0154521074915455,
        -0.016832444605689657, -0.018618915120117006, -0.018760218489373898, -0.01934270325946879,
        -0.019488049955432584, -0.019765249227463824, -0.0206791271808474, -0.02178548323137845,
        -0.022572223061386884, -0.024115903713732928, -0.024275014846415817, -0.02500975830292993,
        -0.025375714345741495, -0.026866214545443433, -0.029468877487032513, -0.03021405313622627,
        -0.03057617559860784, -0.031116322078971602, -0.03171936248193976, -0.03211231301255437,
        -0.0324240720641027, -0.032782765967457966, -0.03425059495071361, -0.03633709430785735,
        -0.038714894661790886, -0.03876518728309634, -0.039082253391589286, -0.03979668659960617,
    ]
    PIN_D_NPAR = [6, 4, 1, 1, 2, 2, 4, 3, 5, 6, 3, 3, 5, 4, 5, 5, 5, 5, 4, 5, 6, 8, 7, 5, 8, 5, 5, 5, 6, 3, 6, 6, 6, 6, 7, 4, 6, 9, 8, 5, 7, 9, 7, 7, 10, 8, 7, 8, 6, 8, 7, 6, 8, 10, 9, 5, 9, 9, 9, 9]
    PIN_D_LIVE = [
        -0.06543881715652791, -0.06135230840369846, -0.059373231315702095, -0.05696407586331669,
        -0.05459929918574495, -0.053872570345624156, -0.051918387390658745, -0.051855255096153664,
        -0.04612657036836043, -0.04460900795433882, -0.043489268748281504, -0.04315188689350818,
    ]

    # the trajectory vectors are asserted on Julia 1.10 only (header); the
    # replay identity and the charge, order, and weld asserts run on every version
    pin3_assert_trajectory = VERSION < v"1.11"

    @testset "fixture S: surface-aware step loop through the empty plateau (seed 424271)" begin
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
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
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
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
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        df, seen, live = pin3_driver_run(424273, pin3_counts_d, 6.0, 40, 60)
        @test nrow(df) == 60
        @test names(df) == ["iter", "emax", "num_particles", "log_compression", "n_live", "n_obs"]  # compression keyword: the ledger gained the n_live column
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
