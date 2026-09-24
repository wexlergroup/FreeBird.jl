# Absolute seeded pins for the atomistic ideal-gas-referenced grand-canonical
# step and kernel, captured at dev 1c7d4494 (the shipped SHA before the
# bounded-support and surface rounds) and reproduced digit-identically across
# two separate Julia processes before recording. The shipped end-to-end
# testsets compare two same-process runs, so no absolute cross-change pin
# existed for this route; stream-neutral changes cite these pins as their
# bit-identity gate.
# Re-recorded for #287 on dev 1130c74e on the per-trial Metropolis uniform stream (every
# insertion and deletion proposal draws its uniform with the proposal, whatever the
# ceiling outcome and the ratio) and again reproduced digit-identically across two
# separate Julia 1.10.4 processes: the particle-count sequences, the emax and live-set
# energies and the kernel counters below are the new stream's, while the fixtures, seeds,
# tolerances, iteration sequences and compression charges (the plateau structure of both
# step-loop fixtures) are unchanged.
#
# Scope: the pins deliberately EXCLUDE the driver's reference-law
# initialization. That path draws its particle counts through Distributions'
# Poisson sampler, whose draw stream is version-dependent (a first capture
# through the initializer diverged from draw one on every CI leg at
# z0V = 6.0 while z0V = 12 survived: a sampler-algorithm boundary), so the
# fixtures build their live sets deterministically (fixed counts, positions
# from raw uniforms) and pin the step-loop stream, which consumes only
# Random-stdlib draws — the family the lattice trajectory pins prove stable
# across CI legs. The initialization law itself is gated by the same-process
# statistical and stream-identity testsets.
#
# Cross-architecture policy (the x64-drift rules): step-loop emax and live-set
# energies are smooth Lennard-Jones accumulations, pinned elementwise at rtol
# 1e-12 (a real stream change moves them by orders of magnitude);
# num_particles sequences are exact integers; compression charges are compared
# against same-process logs of integer ratios (architecture-exact); the
# kernel counter pins run on a zero-interaction fixture where every recorded
# energy is exactly 0.0 eV, so exact equality is defensible on every CI leg.
# Counter asserts are per key and by name, never an exact key set.
#
# Fixture A (seed 424261): K = 16 with two empty walkers, so the descent
# traverses the exact E = 0 plateau (five eviction charges, log(15/16) down
# to log(11/12)). Fixture B (seed 424262): K = 12, denser counts 1..12,
# two eviction charges. Kernel fixture (seed 424253): 500 steps at
# epsilon = 0 under a 1 eV ceiling.
@testset "atomistic igref driver pins (absolute, captured at 1c7d4494)" begin
    using Random

    pin_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
    pin_pbc = (true, true, true)
    pin_seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], pin_box, pin_pbc))
    pin_mkempty() = FastSystem(cell_vectors(pin_seed_at), periodicity(pin_seed_at),
                               empty(position(pin_seed_at, :)), empty(species(pin_seed_at, :)),
                               empty(mass(pin_seed_at, :)))
    pin_V = 1728.0

    function pin_liveset(counts, lj)
        walkers = AtomWalker{1}[]
        for n in counts
            w = AtomWalker{1}(pin_mkempty())
            for _ in 1:n
                pos = SVector(rand() * 12.0, rand() * 12.0, rand() * 12.0)u"Å"
                FreeBird.AbstractWalkers.insert_particle!(w, pos, :Ar)
            end
            push!(walkers, w)
        end
        return GenericAtomWalkers(walkers, lj)
    end

    function pin_step_run(seed::Int, counts, z0V::Float64, mc_steps::Int, n_steps::Int)
        Random.seed!(seed)
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
        ls = pin_liveset(counts, lj)
        params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
            reference_activity=(z0V / pin_V)u"Å^-3", species=:Ar,
            allowed_fail_count=100_000, compression=:mean)  # compression keyword: fixture on the historical mean convention (compression=:mean); the geometric default is covered by test-compression-convention.jl
        iters = Int[]
        emaxs = Float64[]
        npars = Int[]
        logts = Float64[]
        for k in 1:n_steps
            iter, emax, n_par, ls, params, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(
                ls, params, MCAtomGrandCanonicalMoves(); ns_iteration=k, z0V=z0V)
            if !(iter isa Missing)
                push!(iters, iter)
                push!(emaxs, ustrip(u"eV", emax))
                push!(npars, n_par)
                push!(logts, log_t)
            end
        end
        live = sort([ustrip(u"eV", w.energy) for w in ls.walkers])
        return iters, emaxs, npars, logts, live
    end

    function pin_kernel_run(seed::Int)
        Random.seed!(seed)
        w = AtomWalker{1}(pin_mkempty())
        for pos in ([2.0, 2.0, 2.0], [7.0, 7.0, 7.0], [2.0, 7.0, 2.0])
            FreeBird.AbstractWalkers.insert_particle!(w, SVector{3}(pos)u"Å", :Ar)
        end
        w.energy = 0.0u"eV"
        lj0 = LJParameters(epsilon=0.0)
        accept, rate, w, stats = MC_grand_canonical_walk!(500, w, lj0, 1.0u"eV";
            z0V=4.0, species=:Ar, p_move=0.4, p_insert=0.3, step_size=0.8)
        return accept, rate, w.list_num_par[1], stats
    end

    pin_counts_a = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8]
    pin_counts_b = collect(1:12)
    PIN_A_LC_NUM = [16, 16, 15, 14, 13, 12, 11, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16]
    PIN_A_LC_DEN = [17, 17, 16, 15, 14, 13, 12, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17]
    PIN_B_LC_NUM = [12, 12, 12, 12, 12, 11, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
    PIN_B_LC_DEN = [13, 13, 13, 13, 13, 12, 11, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
    PIN_A_EMAX = [
        0.49650308987844566, 0.005567784408978168, 0.0, 0.0,
        0.0, 0.0, 0.0, -0.0001073907988809283,
        -0.00011292767604743266, -0.00020951367630830147, -0.0004954817565344773, -0.0004981137277853385,
        -0.003353093382494646, -0.003619870431342939, -0.004067782910582651, -0.007674231541454843,
        -0.008321397827873075, -0.008615857914117453, -0.008869494002397662, -0.010292890318992425,
        -0.010566885349045572, -0.01077104864357166, -0.011341091625419764, -0.013050278988552816,
        -0.013829662062892805, -0.01446436947367205, -0.014610876566821779, -0.01493769299694203,
        -0.015013997497896149, -0.015918537969139956, -0.017772734750851036, -0.018281583100231084,
        -0.019964122831909355, -0.020219709978764305, -0.020420142085244165, -0.02135002281698043,
        -0.02163454014394159, -0.02166104598262431, -0.02202457017729143, -0.022273956345478324,
        -0.023056425371982596, -0.02325548626049523, -0.02328767122254665, -0.02334588135037035,
        -0.023650510349953194, -0.024348228258367223, -0.02499064848306374, -0.02524699026048351,
        -0.02528623928844231, -0.026505577289071238, -0.027987903006333067, -0.028902932874718422,
        -0.029127273943559708, -0.029612342920634596, -0.03017315285792914, -0.030475962575221633,
        -0.031459120536279195, -0.031476398902901465, -0.03161295552706551, -0.03202937631334223,
        -0.033142190509347906, -0.033687810841705225, -0.0343143314198581, -0.03457167022687381,
        -0.03532173369632021, -0.035416112581165476, -0.03572002314051175, -0.03592464350205273,
        -0.03644510939802923, -0.03653295162293566, -0.03685230130132019, -0.037129424340572224,
        -0.037235943983547695, -0.0390543602118365, -0.03925504496108646, -0.04038091000941275,
        -0.04167259815676241, -0.041687998010135625, -0.04224053922070141, -0.042258970370981995,
        -0.043566773450868616, -0.045145192607996515, -0.04634878500943443, -0.04655358411525691,
        -0.04688769289962708, -0.04908479962672949, -0.04922725347214105, -0.04930398554107008,
        -0.04957932222125763, -0.049730642576068136, -0.04974725079522224, -0.049929083754285494,
        -0.05065896224019509, -0.05096721358332101, -0.05212537173587034, -0.053277396690000596,
        -0.05397186897727904, -0.054644119373167266, -0.05468003286582801, -0.056129687527756725,
        -0.056444042874230806, -0.05668808417527203, -0.05721670304070103, -0.05765946600176303,
        -0.05766026050898085, -0.05806615310940526, -0.0584308663501795, -0.058964366751293754,
        -0.059497408907153757, -0.06094893140592788, -0.06105025592084982, -0.062042613944856224,
        -0.06327297007732842, -0.06351646019037517, -0.0636753921751504, -0.06432652578772463,
        -0.06461574119994935, -0.06472433862228003, -0.06489517319723617, -0.0650087396538078,
    ]
    PIN_A_NPAR = [6, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 3, 4, 3, 3, 3, 8, 5, 4, 6, 4, 6, 5, 9, 7, 7, 4, 3, 7, 7, 5, 8, 7, 7, 5, 9, 9, 8, 5, 7, 7, 8, 6, 5, 5, 7, 7, 7, 8, 8, 7, 8, 8, 9, 8, 7, 9, 7, 10, 5, 7, 9, 6, 8, 8, 7, 9, 9, 8, 7, 9, 7, 8, 9, 10, 10, 9, 8, 11, 7, 7, 7, 9, 11, 9, 10, 9, 9, 9, 11, 9, 8, 11, 8, 6, 9, 9, 8, 8, 12, 9, 9, 10, 9, 13, 10, 11, 10, 10, 11, 8, 12, 9, 12, 11, 10, 9, 12, 11, 12]
    PIN_A_LIVE = [
        -0.08112255104286893, -0.07877062044616945, -0.07674592346866523, -0.07518307918913196,
        -0.07427497348311196, -0.07190176684906507, -0.07179873062942281, -0.07178282179252884,
        -0.07106151152730703, -0.06763657679471814, -0.06723441801879092, -0.06694044368612415,
        -0.06632359620367474, -0.06577678782436684, -0.06557275681864219, -0.06534281338973796,
    ]
    PIN_B_EMAX = [
        2994.8883255499827, 2.5786639977428, 0.3100445290504687, 0.021857781611742764,
        0.015175990895690434, 0.0, 0.0, -0.00017487545907056594,
        -0.0010473248913946223, -0.01024700452087205, -0.012222907546631514, -0.012363430379531987,
        -0.015501215227832115, -0.015662660229791712, -0.021137837313083434, -0.02206647153762751,
        -0.023792074718146324, -0.028541049949380073, -0.02891614831700851, -0.03510556169596466,
        -0.036477755389023, -0.03952228883665778, -0.041711993790707405, -0.04681702596827132,
        -0.04696328800019199, -0.04718475095209054, -0.048270289930232545, -0.04978614208038387,
        -0.05278060049603577, -0.05476921527080571, -0.055139232096943915, -0.055315954743368256,
        -0.055578263062586036, -0.056831792399656716, -0.05720048885997932, -0.05720169667797126,
        -0.05748955465472688, -0.06034810684283498, -0.062038883950396806, -0.06218916908019923,
        -0.06496526280987099, -0.06528882905398879, -0.06589542828704793, -0.06728299748033076,
        -0.07283235963424427, -0.07389591363920636, -0.07456270243507475, -0.07518714022159577,
        -0.07950896421474833, -0.07978705749216719, -0.0812923689842362, -0.082715780392053,
        -0.08379163866301173, -0.08388314815881327, -0.08434903971762875, -0.08594399291206106,
        -0.08764407285442599, -0.08774499576497562, -0.08815723451399683, -0.09030491895390266,
        -0.0909909859240737, -0.09718732538932663, -0.09780455017343606, -0.09839874188625468,
        -0.10129261177173955, -0.10281856774618593, -0.10340950065156801, -0.1045558974195283,
        -0.10463036721870464, -0.10466197267219512, -0.10490901412013537, -0.10535945031781176,
        -0.10583230030922286, -0.10621378309423041, -0.1098979269541191, -0.11101521811252792,
        -0.11310378377885887, -0.11388934359367485, -0.11466650390257435, -0.11546299755204062,
    ]
    PIN_B_NPAR = [11, 8, 8, 6, 12, 1, 2, 3, 4, 7, 5, 8, 7, 7, 12, 10, 7, 8, 8, 10, 10, 10, 11, 14, 12, 10, 11, 11, 13, 12, 10, 11, 10, 12, 11, 10, 10, 12, 12, 13, 12, 11, 13, 15, 14, 9, 12, 13, 13, 13, 14, 14, 14, 14, 12, 14, 14, 15, 15, 12, 13, 15, 14, 15, 15, 15, 12, 13, 15, 14, 17, 13, 16, 14, 16, 15, 15, 15, 16, 15]
    PIN_B_LIVE = [
        -0.15510071803977227, -0.13514983972078334, -0.13210417210173953, -0.13025236640425591,
        -0.12710135927473817, -0.12522846750054184, -0.12519056134598786, -0.12490534700342025,
        -0.12481143883632548, -0.12144202858904363, -0.11965297535223665, -0.11616122274836409,
    ]
    PIN_K_STATS = (move_attempted = 185, move_accepted = 185, insert_attempted = 168,
                   insert_accepted = 120, insert_biased_attempted = 0,
                   insert_biased_accepted = 0, delete_attempted = 139,
                   delete_accepted = 120)

    @testset "fixture A: step-loop descent through the E = 0 plateau (seed 424261)" begin
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        iters, emaxs, npars, logts, live = pin_step_run(424261, pin_counts_a, 6.0, 60, 120)
        @test iters == collect(1:120)
        @test issorted(emaxs, rev=true)
        @test npars == PIN_A_NPAR
        @test logts == log.(PIN_A_LC_NUM ./ PIN_A_LC_DEN)
        @test all(isapprox.(emaxs, PIN_A_EMAX; rtol=1e-12, atol=0.0))
        @test all(isapprox.(live, PIN_A_LIVE; rtol=1e-12, atol=0.0))
    end

    @testset "fixture B: denser step-loop, brief plateau (seed 424262)" begin
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        iters, emaxs, npars, logts, live = pin_step_run(424262, pin_counts_b, 12.0, 40, 80)
        @test iters == collect(1:80)
        @test issorted(emaxs, rev=true)
        @test npars == PIN_B_NPAR
        @test logts == log.(PIN_B_LC_NUM ./ PIN_B_LC_DEN)
        @test all(isapprox.(emaxs, PIN_B_EMAX; rtol=1e-12, atol=0.0))
        @test all(isapprox.(live, PIN_B_LIVE; rtol=1e-12, atol=0.0))
    end

    @testset "kernel counters on the zero-interaction fixture (seed 424253)" begin
        # per-trial Metropolis uniform: re-recorded on the new stream (#287 on dev 1130c74e); fixture, seeds and tolerances unchanged
        accept, rate, n_final, stats = pin_kernel_run(424253)
        @test accept === true
        @test rate == 425 / 500
        @test n_final == 3
        for (k, v) in pairs(PIN_K_STATS)
            @test getproperty(stats, k) == v
        end
    end
end
