# Absolute seeded pins for the atomistic ideal-gas-referenced grand-canonical
# step and kernel, captured at dev 1c7d4494 (the shipped SHA before the
# bounded-support and surface rounds) and reproduced digit-identically across
# two separate Julia processes before recording. The shipped end-to-end
# testsets compare two same-process runs, so no absolute cross-change pin
# existed for this route; stream-neutral changes cite these pins as their
# bit-identity gate.
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
            allowed_fail_count=100_000)
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
        -0.00011292767604743266, -0.00020951367630830147, -0.0008239444072923945, -0.0009805158933565044,
        -0.0015397106957732306, -0.002989816253236788, -0.003353093382494646, -0.006613094585540791,
        -0.0067153573457302915, -0.00714468932937065, -0.007679330050948285, -0.008063262509870175,
        -0.008470218760268954, -0.009919079870469902, -0.01077104864357166, -0.011238323672012163,
        -0.012759464584894964, -0.014610876566821779, -0.015101172692519793, -0.015918537969139956,
        -0.01682128277716249, -0.016914667608249133, -0.017772734750851036, -0.018502451515889733,
        -0.01855795113253103, -0.020416578292040015, -0.02166104598262431, -0.02186733542220379,
        -0.022548321323137867, -0.02268780967025878, -0.02308152338599702, -0.023233069983609626,
        -0.023365278240072864, -0.023373803686156203, -0.023383701360065803, -0.023811600490849908,
        -0.024514675928337223, -0.02462332825424484, -0.02541383387245483, -0.025885173425747938,
        -0.02623504028045993, -0.0264979433966434, -0.02679504293639932, -0.02752305688172925,
        -0.028988678482161762, -0.029668629317186855, -0.029890324681572374, -0.03037292210685511,
        -0.03198113310138711, -0.032380163888193224, -0.03255871294899312, -0.032968235701439925,
        -0.032997717356232706, -0.033101585522182876, -0.03357304468058118, -0.03431780282923257,
        -0.035189526107308924, -0.036274045873984716, -0.036956719888022634, -0.03721872739087561,
        -0.03797619945797363, -0.03933189235986552, -0.03989052741524547, -0.040181033217279304,
        -0.04089523694188833, -0.0414046508333025, -0.041910551078641556, -0.04196269332817865,
        -0.04250006320498496, -0.04255831377183063, -0.042638674284410115, -0.04289385290378052,
        -0.04330866455818437, -0.04396193365349088, -0.04498262635237586, -0.04511044710485016,
        -0.04576681639184778, -0.04619034410132479, -0.046549818014362546, -0.04667064243790245,
        -0.047168077399936115, -0.04756524728975553, -0.05046584065965944, -0.05110468871686705,
        -0.05113774207550152, -0.051853168551060196, -0.05243704402998366, -0.053015801846527655,
        -0.053098611099561446, -0.05362580987820714, -0.05384980781518607, -0.054357147349725576,
        -0.05574277702802663, -0.05654864560128079, -0.056641132542510575, -0.05666317770622394,
        -0.056698811202237825, -0.056995386106409994, -0.05709415610764919, -0.05718716019707145,
        -0.05758342515432109, -0.05889880077227439, -0.05943441830518491, -0.05965605862959167,
        -0.06017244074537433, -0.06033673204260611, -0.060357743462356415, -0.06047541259886642,
        -0.061292457420684406, -0.06178768802860045, -0.06217256714590098, -0.06227598114412991,
    ]
    PIN_A_NPAR = [6, 5, 0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 4, 5, 5, 5, 7, 5, 4, 6, 6, 7, 8, 4, 7, 7, 7, 6, 5, 4, 7, 8, 8, 7, 9, 8, 8, 8, 8, 7, 8, 5, 6, 6, 8, 7, 8, 6, 6, 7, 9, 9, 5, 8, 7, 8, 9, 9, 9, 10, 7, 6, 9, 10, 7, 9, 11, 8, 9, 9, 9, 9, 10, 10, 9, 10, 10, 11, 8, 8, 8, 8, 11, 9, 10, 11, 8, 10, 6, 12, 11, 10, 10, 12, 12, 10, 10, 9, 13, 9, 11, 11, 9, 13, 10, 11, 10, 13, 10, 11, 12, 9, 10, 7, 11, 10, 12, 10]
    PIN_A_LIVE = [
        -0.09054368869280663, -0.08154836655860949, -0.08029997698197638, -0.07932163124190357,
        -0.07647101166652966, -0.07414024477961706, -0.07197632266225525, -0.0718033058069537,
        -0.06735075149873104, -0.06630696075682387, -0.06537782030975144, -0.06481040922359776,
        -0.06476584940293108, -0.0642550289296039, -0.06395595518186165, -0.06361822740367877,
    ]
    PIN_B_EMAX = [
        2994.8883255499827, 2.5786639977428, 0.021857781611742764, 0.015175990895690434,
        0.010466702437576894, 0.0, 0.0, -0.00017487545907056594,
        -0.0010473248913946223, -0.0013015014148710763, -0.006080130314878283, -0.012222907546631514,
        -0.013907433160176714, -0.015662660229791712, -0.0197045361615375, -0.020227452313412064,
        -0.02206647153762751, -0.025531565917946974, -0.027050850745844762, -0.02855588273989433,
        -0.028690240502911024, -0.030729841846923154, -0.03135437380029623, -0.03654248956410261,
        -0.038119721925529, -0.040839973188740046, -0.04238728993424989, -0.04310937864021593,
        -0.04865735378455641, -0.048974415318777005, -0.049608057997187954, -0.05059789513763672,
        -0.05542170045583214, -0.05737827840565489, -0.05763174410110945, -0.059334237443026926,
        -0.060571515378430925, -0.061071620458641986, -0.06141775568292222, -0.062124235607597444,
        -0.06582601815935438, -0.06662380944067647, -0.06770403589364778, -0.06916599340910214,
        -0.06985532997001528, -0.07065855880587775, -0.07264138273093595, -0.07374240093506991,
        -0.07389591363920636, -0.07414890724572859, -0.07523405635410964, -0.07661264565963569,
        -0.07676425158393256, -0.07947957566465773, -0.08096480049231033, -0.08302045083962474,
        -0.08329833535651876, -0.08443055225384373, -0.08464043929604013, -0.08795732383242695,
        -0.08817614171131019, -0.09154960665147996, -0.09217625435213611, -0.09219713969106322,
        -0.09232788252572481, -0.09359550720145389, -0.09735988456989661, -0.09872732489426503,
        -0.10187954523310573, -0.10237607482008534, -0.10531050081049796, -0.10583082651691188,
        -0.10675036971327378, -0.10697625754129211, -0.11041761737417012, -0.11131593516741708,
        -0.11179657708950429, -0.11303993835732476, -0.11344219554342466, -0.11448433036214237,
    ]
    PIN_B_NPAR = [11, 8, 6, 12, 10, 1, 2, 3, 4, 7, 6, 5, 7, 7, 7, 9, 10, 11, 7, 9, 10, 8, 8, 9, 11, 9, 13, 10, 10, 10, 14, 11, 12, 9, 11, 14, 10, 11, 11, 12, 11, 11, 13, 13, 13, 11, 13, 13, 9, 12, 12, 13, 16, 14, 14, 13, 15, 14, 15, 13, 14, 13, 14, 11, 17, 13, 14, 13, 16, 14, 15, 15, 16, 14, 14, 15, 13, 17, 14, 14]
    PIN_B_LIVE = [
        -0.1298580210467591, -0.12823039993363036, -0.12566080955215025, -0.12241058656469743,
        -0.12166589708845912, -0.1201598350549118, -0.1192737383086962, -0.11883694052493093,
        -0.11797676514706185, -0.11708551033046782, -0.11603884659794568, -0.11590144088572349,
    ]
    PIN_K_STATS = (move_attempted = 181, move_accepted = 181, insert_attempted = 155,
                   insert_accepted = 132, insert_biased_attempted = 0,
                   insert_biased_accepted = 0, delete_attempted = 159,
                   delete_accepted = 135)

    @testset "fixture A: step-loop descent through the E = 0 plateau (seed 424261)" begin
        iters, emaxs, npars, logts, live = pin_step_run(424261, pin_counts_a, 6.0, 60, 120)
        @test iters == collect(1:120)
        @test issorted(emaxs, rev=true)
        @test npars == PIN_A_NPAR
        @test logts == log.(PIN_A_LC_NUM ./ PIN_A_LC_DEN)
        @test all(isapprox.(emaxs, PIN_A_EMAX; rtol=1e-12, atol=0.0))
        @test all(isapprox.(live, PIN_A_LIVE; rtol=1e-12, atol=0.0))
    end

    @testset "fixture B: denser step-loop, brief plateau (seed 424262)" begin
        iters, emaxs, npars, logts, live = pin_step_run(424262, pin_counts_b, 12.0, 40, 80)
        @test iters == collect(1:80)
        @test issorted(emaxs, rev=true)
        @test npars == PIN_B_NPAR
        @test logts == log.(PIN_B_LC_NUM ./ PIN_B_LC_DEN)
        @test all(isapprox.(emaxs, PIN_B_EMAX; rtol=1e-12, atol=0.0))
        @test all(isapprox.(live, PIN_B_LIVE; rtol=1e-12, atol=0.0))
    end

    @testset "kernel counters on the zero-interaction fixture (seed 424253)" begin
        accept, rate, n_final, stats = pin_kernel_run(424253)
        @test accept === true
        @test rate == 448 / 500
        @test n_final == 0
        for (k, v) in pairs(PIN_K_STATS)
            @test getproperty(stats, k) == v
        end
    end
end
