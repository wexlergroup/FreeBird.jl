# The serial atomistic nested-sampling steps keep the clone of a decorrelation walk that
# accepted no move (`_unmoved_clone_below`): the canonical random-walk steps (3D and 2D),
# the surface step, the mixed-moves step, the Galilean step and the two
# ideal-gas-referenced grand-canonical steps, in the ordinary cull and in the plateau
# refill. Every rejected proposal of the atomistic kernels restores the previous state
# exactly, so such a clone is its parent unchanged: a survivor below the ceiling, and a
# sample of the constrained prior like any walked clone. Retrying from another parent
# instead would copy walkers whose walks rarely move too rarely. Below, every walk
# accepts no move: with `mc_steps = 0`, or (the canonical random-walk steps, whose rate
# bookkeeping needs a step) with one displacement of up to 10^6 Å across a slit along x,
# which lands outside the field's accessible region and is rejected.
#
# A kept clone ties its parent exactly. A block of such copies is a tie of copying and not an
# atom of the energy law, so its evictions are charged as ordinary culls at the block's
# starting live count (`_copy_tie_block`), while the empty configuration and tied walkers
# whose configurations differ keep the atom charge ln((n - 1)/n).
#
# Calibration ledger (the gates ship at >= 3x the largest deviation over the
# calibration blocks, per statistic):
# - Reference-fugacity descent of the ideal gas in the planar slit U = 0.005 (z - 1) eV/Å
#   on [1, 9] Å (K = 256, z0 V_acc = 4, four-step walks, a block of 32 seeded runs,
#   ln Ξ and <N> at the reference activity and 300 K against z I, exact <N> = 2.035): the
#   largest |block mean deviation| over twelve blocks (seeds 1 to 384) is 0.027 in ln Ξ
#   and 0.043 in <N>; the gates ship at 0.09 and 0.13. Over the same twelve blocks the
#   previous rule (discard the clone and retry from another parent) deviates in <N> by
#   +0.196 to +0.236, every block outside the gate (+0.206 on the shipped block); its
#   ln Ξ deviations, -0.022 to -0.064, stay inside that gate, so <N> is the
#   discriminating statistic. The shipped block is seeds 1 to 32.
# - One-move canonical descent of one particle in U = 0.05 (z - 1) eV/Å on [1, 9] Å
#   (K = 64), charged ln X against the closed form at depth 15, a block of 32 seeded runs:
#   the largest |block mean deviation| over twelve blocks (seeds 1 to 384) is 1.13; the gate
#   ships at 3.5. Charging its copy blocks with the atom factor instead deviates by -11.6
#   to -15.0, every block outside the gate. The shipped block is seeds 1 to 32.
# - Every other testset is an exact contract and needs no calibration.
@testset "Atomistic NS steps keep an unmoved clone" begin
    using Random

    # >>> fixtures (the calibration script extracts this block verbatim)
    ac_kb = 8.617333262e-5
    ac_mass = 39.948u"u"
    ac_lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(ac_mass, T * u"K"))
    ac_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    ac_slit_at = FastSystem(atomic_system([:Ar => [6.0, 6.0, 5.0]u"Å"], ac_box, (true, true, false)))
    # an ideal gas in the linear planar field U = 0.005 (z - 1) eV/Å on [1, 9] Å
    ac_slit = TabulatedPlanarField(3, [1.0, 9.0]u"Å", [0.0, 0.04]u"eV")
    ac_pot = ExternalFieldPotential(IdealGasParameters(), ac_slit)
    ac_mkempty() = FastSystem(cell_vectors(ac_slit_at), periodicity(ac_slit_at),
                              empty(position(ac_slit_at, :)), empty(species(ac_slit_at, :)),
                              empty(mass(ac_slit_at, :)))
    ac_save(tag) = SaveEveryN(df_filename="_unmoved_$(tag).csv",
                              wk_filename="_unmoved_$(tag).traj.extxyz",
                              ls_filename="_unmoved_$(tag).ls.extxyz",
                              n_traj=10^9, n_snap=10^9, n_info=10^9)
    ac_clean(tag) = for f in ["_unmoved_$(tag).csv", "_unmoved_$(tag).traj.extxyz",
                              "_unmoved_$(tag).ls.extxyz"]
        rm(f, force=true)
    end
    # I(T) = ∫ e^{-βU} dr over the slit (area 144 Å^2, width 8 Å)
    ac_I(T) = (β = 1 / (ac_kb * T); 144.0 * (1 - exp(-β * 0.04)) / (β * 0.005))
    # one reference-fugacity descent; deviations of ln Ξ and <N> at the reference
    # activity and temperature T from the closed forms ln Ξ = <N> = z I
    function ac_igref_dev(seed; K=256, mc_steps=4, z0V_acc=4.0, n_steps=20_000, T=300.0)
        Random.seed!(seed)
        ls = GenericAtomWalkers([AtomWalker{1}(ac_mkempty()) for _ in 1:K], ac_pot)
        Vacc = accessible_volume(ac_slit, ac_slit_at)
        z0 = z0V_acc / Vacc
        params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps, reference_activity=z0,
                                              species=:Ar, allowed_fail_count=1000)
        # the run ends on the empty configuration (the stall stop): silence its warnings
        df, lso, _ = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            ideal_gas_referenced_nested_sampling(ls, params, n_steps, MCAtomGrandCanonicalMoves(), ac_save("o"))
        end
        ac_clean("o")
        z = ustrip(u"Å^-3", z0)
        mu = (ac_kb * T * (log(z) + 3 * log(ac_lam(T)))) * u"eV"
        st = gc_thermodynamic_stats_ideal_ref(df, Vacc, ac_mass, z0, [mu], [T]u"K";
            live_emax=[ustrip(u"eV", w.energy) for w in lso.walkers],
            live_numbers=[w.list_num_par[1] for w in lso.walkers])
        return st.logXi[1, 1] - z * ac_I(T), st.mean_N[1, 1] - z * ac_I(T)
    end
    # one particle in U = 0.05 (z - 1) eV/Å on [1, 9] Å: the prior mass below E is E/(0.4 eV).
    # One-move walks copy heavily; the deviation is the charged ln X at the first row past
    # depth 15 minus the closed form
    function ac_copy_dev(seed; K=64, n_steps=1600, depth=15.0, Ec=0.4)
        Random.seed!(seed)
        pot = ExternalFieldPotential(IdealGasParameters(), TabulatedPlanarField(3, [1.0, 9.0]u"Å", [0.0, Ec]u"eV"))
        walkers = map(1:K) do _
            w = AtomWalker{1}(ac_mkempty())
            FreeBird.AbstractWalkers.insert_particle!(w, SVector(12 * rand(), 12 * rand(), 1 + 8 * rand())u"Å", :Ar)
            AtomWalker{1}(w.configuration)
        end
        params = NestedSamplingParameters(mc_steps=1, step_size=2.0, initial_step_size=2.0, allowed_fail_count=10^7)
        df, _, _ = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            nested_sampling(GenericAtomWalkers(walkers, pot), params, n_steps, MCRandomWalkClone(), ac_save("cp"))
        end
        ac_clean("cp")
        i = findlast(>=(Ec * exp(-depth)), df.emax)
        return sum(df.log_compression[1:i]) - log(df.emax[i] / Ec)
    end
    # <<< fixtures

    # a walker of the slit with one particle per height in zs (in-plane positions fixed)
    function ac_field_walker(zs)
        w = AtomWalker{1}(ac_mkempty())
        for (k, z) in enumerate(zs)
            FreeBird.AbstractWalkers.insert_particle!(w, SVector(1.0 + k, 6.0, z)u"Å", :Ar)
        end
        w.energy = interacting_energy(w.configuration, ac_pot, w.list_num_par, w.frozen)
        return w
    end
    ac_field_ls(heights) = GenericAtomWalkers([ac_field_walker(zs) for zs in heights], ac_pot)
    # the same slit along x, for the canonical random-walk steps (3D and in-plane walks)
    ac_x_at = FastSystem(atomic_system([:Ar => [5.0, 6.0, 5.0]u"Å"], ac_box, (false, true, true)))
    ac_xpot = ExternalFieldPotential(IdealGasParameters(), TabulatedPlanarField(1, [1.0, 9.0]u"Å", [0.0, 0.04]u"eV"))
    function ac_x_walker(x)
        w = AtomWalker{1}(FastSystem(atomic_system([:Ar => [x, 6.0, 5.0]u"Å"], ac_box, (false, true, true))))
        w.energy = interacting_energy(w.configuration, ac_xpot, w.list_num_par, w.frozen)
        return w
    end
    ac_x_ls(xs) = GenericAtomWalkers([ac_x_walker(x) for x in xs], ac_xpot)
    # one displacement of up to 10^6 Å per walk: it leaves the accessible region
    bigp() = NestedSamplingParameters(mc_steps=1, step_size=1e6, initial_step_size=1e6,
                                      step_size_up=1e7, allowed_fail_count=1000)
    # Lennard-Jones pairs in a periodic cubic cell (the Galilean kernel needs one)
    ac_L = 12.0
    ac_cube = [[ac_L, 0.0, 0.0], [0.0, ac_L, 0.0], [0.0, 0.0, ac_L]]u"Å"
    ac_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.0, shift=true)
    ac_pair(r) = AtomWalker(FastSystem(periodic_system(
        [:Ar => [0.3, 0.5, 0.5], :Ar => [0.3 + r / ac_L, 0.5, 0.5]], ac_cube, fractional=true)))
    # one adsorbate above a frozen four-atom surface
    sf_box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]]u"Å"
    sf_surf_sys = FastSystem(atomic_system(
        [:H => [2.5, 2.5, 2.0]u"Å", :H => [7.5, 2.5, 2.0]u"Å",
         :H => [2.5, 7.5, 2.0]u"Å", :H => [7.5, 7.5, 2.0]u"Å"], sf_box, (true, true, false)))
    sf_cps = CompositeParameterSets(2, [LJParameters(epsilon=0.001, sigma=2.5, cutoff=1.8, shift=true),
                                        LJParameters(epsilon=0.003, sigma=2.5, cutoff=1.8, shift=true),
                                        LJParameters(epsilon=0.01, sigma=2.5, cutoff=1.8, shift=true)])
    sf_walker(z) = AtomWalker{1}(FastSystem(atomic_system([:H => [3.0, 3.0, z]u"Å"], sf_box, (true, true, false))))
    sf_ls(zs) = LJSurfaceWalkers([sf_walker(z) for z in zs], sf_cps,
                                 AtomWalker(deepcopy(sf_surf_sys); freeze_species=[:H]); assign_energy=true)

    snap(w) = (copy(position(w.configuration, :)), w.list_num_par[1], w.energy)
    is_copy_of(w, s) = snap(w) == s
    nsp(mc_steps) = NestedSamplingParameters(mc_steps=mc_steps, step_size=0.5, initial_step_size=0.5,
                                             allowed_fail_count=1000)
    igp(mc_steps; kw...) = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
        reference_activity=(4.0 / 1152.0)u"Å^-3", species=:Ar, allowed_fail_count=1000; kw...)
    igp_sf(mc_steps) = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
        reference_activity=(3.0 / 1500.0)u"Å^-3", species=:H, allowed_fail_count=1000)

    # the cases: (name, live-set builder from six distinct keys, parameters, routine, Ω key)
    canon_heights = [[1.5 + 1.1k] for k in 1:6]   # one particle per walker, along z
    cases = [
        ("random walk, 3D", () -> ac_x_ls([1.5 + 1.1k for k in 1:6]), bigp, MCRandomWalkClone()),
        ("random walk, 2D", () -> ac_x_ls([1.5 + 1.1k for k in 1:6]), bigp, MCRandomWalkClone(dims=[1, 2])),
        ("mixed moves", () -> ac_field_ls(canon_heights), () -> nsp(0), MCMixedMoves(5, 1)),
        ("Galilean walk", () -> LJAtomWalkers([ac_pair(2.80 + 0.03k) for k in 1:6], ac_lj), () -> nsp(0), MCGalileanWalk()),
        ("surface random walk", () -> sf_ls([3.6 + 0.4k for k in 1:6]), () -> nsp(0), MCRandomWalkClone()),
        ("ideal-gas-referenced", () -> ac_field_ls([[1.2 + 0.7j for j in 1:k] for k in 1:6]), () -> igp(0), MCAtomGrandCanonicalMoves()),
        ("ideal-gas-referenced, surface", () -> sf_ls([3.6 + 0.4k for k in 1:6]), () -> igp_sf(0), MCAtomGrandCanonicalMoves()),
    ]
    step!(ls, p, r) = FreeBird.SamplingSchemes.nested_sampling_step!(ls, p, r)

    @testset "_unmoved_clone_below" begin
        f = FreeBird.SamplingSchemes._unmoved_clone_below
        @test f(false, 0.0, 1.0u"eV", 2.0u"eV")
        @test f(false, NaN, 1.0u"eV", 2.0u"eV")          # a walk of zero steps
        @test !f(true, 0.5, 1.0u"eV", 2.0u"eV")          # the walk moved
        @test !f(false, 0.25, 1.0u"eV", 2.0u"eV")        # moved, rejected afterwards
        @test !f(false, 0.0, 2.0u"eV", 2.0u"eV")         # tied with the ceiling
        @test !f(false, 0.0, nextfloat(2.0)u"eV", 2.0u"eV")
    end

    @testset "ordinary cull keeps the unmoved clone: $name" for (name, mk, mkp, routine) in cases
        Random.seed!(7101)
        ls = mk()
        igref = routine isa MCAtomGrandCanonicalMoves
        sort!(ls.walkers, by=w -> w.energy, rev=true)
        ceiling = ls.walkers[1].energy
        @test length(Set(w.energy for w in ls.walkers)) == 6   # six distinct keys
        worst_n = ls.walkers[1].list_num_par[1]
        survivors = [snap(w) for w in ls.walkers[2:end]]
        p = mkp()
        out = step!(ls, p, routine)
        iter, emax = out[1], out[2]
        @test iter isa Int
        @test emax == ceiling
        @test length(ls.walkers) == 6
        kept = ls.walkers[end]
        @test any(s -> is_copy_of(kept, s), survivors)
        @test kept.energy < ceiling
        @test p.fail_count == 1
        if routine isa MCMixedMoves
            @test length(out) == 4                     # this step records no charge
        else
            @test out[end] == -1 / 6                   # an ordinary cull's charge
        end
        igref && @test out[3] == worst_n               # the culled walker's particle count
    end

    @testset "a clone tied with the ceiling is not kept" begin
        # MCRandomWalkMaxE walks a copy of the culled walker itself
        Random.seed!(7102)
        ls = ac_x_ls([1.5 + 1.1k for k in 1:6])
        before = sort([snap(w) for w in ls.walkers])
        p = bigp()
        iter, emax, _, _, log_t = step!(ls, p, MCRandomWalkMaxE())
        @test iter === missing && emax === missing && log_t === missing
        @test sort([snap(w) for w in ls.walkers]) == before
        @test p.fail_count == 1
        # a live set tied in its entirety (the all-tied fallback)
        ls = ac_x_ls(fill(4.0, 4))
        p = bigp()
        iter, emax, _, _, log_t = step!(ls, p, MCRandomWalkClone())
        @test iter === missing && emax === missing
        @test length(ls.walkers) == 4
        ls = GenericAtomWalkers([ac_field_walker([4.0]) for _ in 1:4], ac_pot)
        p = igp(0)
        iter, emax, n, _, _, log_t = step!(ls, p, MCAtomGrandCanonicalMoves())
        @test iter === missing && emax === missing
    end

    # the cases for the refill: the top two walkers are bit-exact duplicates
    refill_cases = [
        ("random walk, 3D", () -> ac_x_ls([7.5, 7.5, (1.5 + 1.1k for k in 1:4)...]), bigp, MCRandomWalkClone()),
        ("random walk, 2D", () -> ac_x_ls([7.5, 7.5, (1.5 + 1.1k for k in 1:4)...]), bigp, MCRandomWalkClone(dims=[1, 2])),
        ("Galilean walk", () -> LJAtomWalkers([ac_pair(2.70), ac_pair(2.70), [ac_pair(2.80 + 0.03k) for k in 1:4]...], ac_lj), () -> nsp(0), MCGalileanWalk()),
        ("surface random walk", () -> sf_ls([3.3, 3.3, (3.6 + 0.4k for k in 1:4)...]), () -> nsp(0), MCRandomWalkClone()),
        ("ideal-gas-referenced", () -> ac_field_ls([[1.2 + 0.7j for j in 1:6], [1.2 + 0.7j for j in 1:6], [[1.2 + 0.7j for j in 1:k] for k in 1:4]...]), () -> igp(0), MCAtomGrandCanonicalMoves()),
        ("ideal-gas-referenced, surface", () -> sf_ls([3.3, 3.3, (3.6 + 0.4k for k in 1:4)...]), () -> igp_sf(0), MCAtomGrandCanonicalMoves()),
    ]
    @testset "plateau refill keeps the unmoved clones: $name" for (name, mk, mkp, routine) in refill_cases
        Random.seed!(7103)
        ls = mk()
        sort!(ls.walkers, by=w -> w.energy, rev=true)
        tie = ls.walkers[1].energy
        @test ls.walkers[2].energy == tie && ls.walkers[3].energy < tie
        survivors = [snap(w) for w in ls.walkers[3:end]]
        p = mkp()
        p.allowed_fail_count = 1   # a kept clone is no failed attempt, so one attempt suffices
        out1 = step!(ls, p, routine)
        @test out1[end] == -1 / 6                      # a block of copies: an ordinary cull's charge
        @test p.plateau_copy_block
        @test length(ls.walkers) == 5
        # the last tied walker: evict it, then refill with copies of the survivors; no
        # refill attempt fails, so nothing is logged
        out2 = @test_logs step!(ls, p, routine)
        @test out2[end] == -1 / 6                      # at the block's starting live count
        @test length(ls.walkers) == 6
        @test p.plateau_refill_target == 0
        @test !p.plateau_copy_block
        @test all(w -> w.energy < tie, ls.walkers)
        @test all(w -> any(s -> is_copy_of(w, s), survivors), ls.walkers)
    end

    @testset "a kept clone and its parent are evicted as a block of exact ties" begin
        # three walkers; the first step culls the worst and keeps a copy of a random
        # survivor. Seeds are tried until the copy is of the middle walker, which it ties
        # exactly, so the case holds on any random stream; the pair is then evicted as a tie
        # block and the refill copies the last survivor
        local ls, mid, p, out
        for seed in 7106:7169
            Random.seed!(seed)
            ls = ac_x_ls([8.0, 5.0, 2.0])
            mid = snap(ls.walkers[2])
            p = bigp()
            out = step!(ls, p, MCRandomWalkClone())
            count(w -> is_copy_of(w, mid), ls.walkers) == 2 && break
        end
        @test out[end] == -1 / 3
        @test count(w -> is_copy_of(w, mid), ls.walkers) == 2
        out = step!(ls, p, MCRandomWalkClone())
        @test out[2] == mid[3] && out[end] == -1 / 3   # copies: ordinary culls from three walkers
        out = step!(ls, p, MCRandomWalkClone())
        @test out[2] == mid[3] && out[end] == -1 / 3
        @test length(ls.walkers) == 3
        @test all(w -> ustrip(u"eV", w.energy) ≈ 0.005, ls.walkers)
    end

    @testset "every step records a row when walks accept no move (canonical driver)" begin
        Random.seed!(7106)
        ls = ac_x_ls([1.3 + 0.8k + 0.01k^2 for k in 1:8])
        df, lso, p = nested_sampling(ls, bigp(), 6, MCRandomWalkClone(), ac_save("drv"))
        ac_clean("drv")
        @test nrow(df) == 6
        @test df.iter == collect(1:6)
        @test issorted(df.emax, rev=true)
        # a kept clone's walk counts towards the driver's step-size reset
        ls = ac_x_ls([1.3 + 0.8k + 0.01k^2 for k in 1:8])
        p0 = bigp()
        p0.allowed_fail_count = 1
        df, lso, p = @test_logs (:warn, r"Reset step size") match_mode=:any nested_sampling(
            ls, p0, 6, MCRandomWalkClone(), ac_save("drv"))
        ac_clean("drv")
        @test p.step_size == p.initial_step_size
    end

    @testset "the rows of a tie block record the live count they cull from" begin
        ls = ac_field_ls([[1.2 + 0.7j for j in 1:6], [1.2 + 0.7j for j in 1:6], [[1.2 + 0.7j for j in 1:k] for k in 1:4]...])
        df, lso, p = ideal_gas_referenced_nested_sampling(ls, igp(0), 3, MCAtomGrandCanonicalMoves(), ac_save("nl"); initialize=false)
        ac_clean("nl")
        @test df.n_live == [6, 5, 6]
        @test df.log_compression == fill(-1 / 6, 3)    # copies, then an ordinary cull
    end

    @testset "_copy_tie_block" begin
        f = FreeBird.SamplingSchemes._copy_tie_block
        w(x) = ac_x_walker(x)
        @test f([w(4.0), w(4.0), w(2.0)], 2)
        @test f([w(4.0), w(4.0), w(4.0)], 3)
        @test !f([w(4.0), w(4.0)], 1)                  # a unique ceiling
        @test !f([w(4.0), w(nextfloat(4.0))], 2)       # one coordinate differs by one ulp
        @test !f([w(4.0), w(4.0), w(nextfloat(4.0))], 3) # a mixed block is an atom
        # two free-particle counts
        @test !f([ac_field_walker([2.0]), ac_field_walker([2.0, 3.0])], 2)
        # two components: the same positions under a different partition are a different state
        two(parts) = AtomWalker{2}(FastSystem(atomic_system(
            [:Ar => [3.0, 3.0, 5.0]u"Å", :Ar => [4.0, 3.0, 5.0]u"Å"], sf_box, (true, true, false)));
            list_num_par=parts, frozen=[false, false])
        @test !f([two([1, 1]), two([2, 0])], 2)
        # the empty configuration is an atom however it arose
        @test !f([AtomWalker{1}(ac_mkempty()), AtomWalker{1}(ac_mkempty())], 2)
        # so is an empty adsorbate layer over a frozen substrate
        subs() = AtomWalker{2}(deepcopy(sf_surf_sys); list_num_par=[0, 4], frozen=[false, true])
        @test !f([subs(), subs()], 2)
        # a substrate holding a free particle in copies is a block of copies
        ads() = AtomWalker{2}(FastSystem(atomic_system(
            [:Ar => [3.0, 3.0, 5.0]u"Å", :H => [2.5, 2.5, 2.0]u"Å"], sf_box, (true, true, false)));
            list_num_par=[1, 1], frozen=[false, true])
        @test f([ads(), ads()], 2)
    end

    # the same six live sets as the keep-rule cases; the two highest keys are made equal by
    # hand on walkers whose configurations differ, as on the vacuum plateau of a truncated potential
    atom_cases = [c for c in cases if !(c[4] isa MCMixedMoves)]
    @testset "tied walkers with different configurations are an atom: $name" for (name, mk, mkp, routine) in atom_cases
        Random.seed!(7112)
        ls = mk()
        sort!(ls.walkers, by=w -> w.energy, rev=true)
        ls.walkers[2].energy = ls.walkers[1].energy
        p = mkp()
        out = step!(ls, p, routine)
        @test out[end] == log(5 / 6)
        @test !p.plateau_copy_block
        out = step!(ls, p, routine)
        @test out[end] == log(4 / 5)
    end

    @testset "a mixed block keeps the atom charge when only copies remain: $name" for (name, mk, mkp, routine) in refill_cases
        Random.seed!(7113)
        ls = mk()
        sort!(ls.walkers, by=w -> w.energy, rev=true)
        # the first walker evicted has a different configuration and the energy of the two
        # bit-exact duplicates behind it (the sorts are stable)
        other = popat!(ls.walkers, 3)
        other.energy = ls.walkers[1].energy
        pushfirst!(ls.walkers, other)
        p = mkp()
        charges = [step!(ls, p, routine)[end] for _ in 1:3]
        @test charges == [log(5 / 6), log(4 / 5), log(3 / 4)]
    end

    @testset "the empty configuration keeps the atom charge" begin
        # at a positive chemical potential every occupied walker lies below Ω = 0
        Random.seed!(7114)
        empties = [AtomWalker{1}(ac_mkempty()) for _ in 1:2]
        ls = GenericAtomWalkers([empties..., (ac_field_walker([1.2 + 0.7j for j in 1:k]) for k in 1:4)...], ac_pot)
        p = igp(0; chemical_potential=0.1u"eV")
        out = step!(ls, p, MCAtomGrandCanonicalMoves())
        @test out[end] == log(5 / 6)
        @test !p.plateau_copy_block
        out = step!(ls, p, MCAtomGrandCanonicalMoves())
        @test out[end] == log(4 / 5)
    end

    @testset "the classification is per block" begin
        # copies at the top, then two tied walkers with different configurations
        Random.seed!(7115)
        ls = ac_x_ls([8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0])
        ls.walkers[3].energy = ls.walkers[4].energy
        p = bigp()
        charges = [step!(ls, p, MCRandomWalkClone())[end] for _ in 1:2]
        @test charges == fill(-1 / 8, 2)
        @test !p.plateau_copy_block
        @test step!(ls, p, MCRandomWalkClone())[end] == log(7 / 8)
    end

    @testset "a block of copies under the mean convention" begin
        Random.seed!(7116)
        ls = ac_x_ls([7.5, 7.5, 4.0, 3.0, 2.0, 1.5])
        p = NestedSamplingParameters(mc_steps=1, step_size=1e6, initial_step_size=1e6,
                                     step_size_up=1e7, allowed_fail_count=1000, compression=:mean)
        charges = [step!(ls, p, MCRandomWalkClone())[end] for _ in 1:2]
        @test charges == fill(log(6 / 7), 2)
    end

    @testset "a kept clone passes through the Galilean burst" begin
        kin(w) = interacting_energy(w.configuration, ac_lj, w.list_num_par, w.frozen)
        routine = MCAtomGrandCanonicalMoves(galilean_steps=4)
        # ordinary cull: the walk accepts nothing and the burst still runs
        Random.seed!(7117)
        ls = LJAtomWalkers([ac_pair(2.80 + 0.03k) for k in 1:6], ac_lj)
        sort!(ls.walkers, by=w -> w.energy, rev=true)
        ceiling = ls.walkers[1].energy
        p = igp(0)
        out = step!(ls, p, routine)
        @test out[1] isa Int && out[end] == -1 / 6
        @test p.move_stats[:galilean_attempted] > 0
        @test p.fail_count == 1
        kept = ls.walkers[end]
        @test kept.energy < ceiling
        @test kept.energy == kin(kept)                 # re-anchored from scratch
        # the plateau refill
        Random.seed!(7118)
        ls = LJAtomWalkers([ac_pair(2.70), ac_pair(2.70), [ac_pair(2.80 + 0.03k) for k in 1:4]...], ac_lj)
        sort!(ls.walkers, by=w -> w.energy, rev=true)
        tie = ls.walkers[1].energy
        survivors = [snap(w) for w in ls.walkers[3:end]]
        p = igp(0)
        step!(ls, p, routine)
        step!(ls, p, routine)
        @test length(ls.walkers) == 6
        @test all(w -> w.energy < tie, ls.walkers)
        @test any(w -> !any(s -> is_copy_of(w, s), survivors), ls.walkers)   # a burst moved a refilled clone
    end

    @testset "the stall stop keeps the row of the step that reached it" begin
        # every walk accepts no move and every clone is kept, so fail_count reaches
        # allowed_fail_count = 3 on a step that recorded a row; that row is in the ledger
        Random.seed!(7107)
        ls = ac_field_ls([[1.2 + 0.7j for j in 1:k] for k in 1:8])
        p = AtomisticIGRefGCNSParameters(mc_steps=0, reference_activity=(4.0 / 1152.0)u"Å^-3",
                                         species=:Ar, allowed_fail_count=3)
        df, lso, p = @test_logs (:warn, r"stop_on_stall") match_mode=:any ideal_gas_referenced_nested_sampling(
            ls, p, 50, MCAtomGrandCanonicalMoves(), ac_save("stall"); initialize=false)
        ac_clean("stall")
        @test p.fail_count == 3
        @test nrow(df) >= 3
        @test df.iter == collect(1:nrow(df))
        @test all(w -> w.iter == nrow(df), lso.walkers)
        @test length(lso.walkers) == 8
        # with the stop off, the run warns, resets and carries on to its last step
        ls = ac_field_ls([[1.2 + 0.7j for j in 1:k] for k in 1:8])
        p = AtomisticIGRefGCNSParameters(mc_steps=0, reference_activity=(4.0 / 1152.0)u"Å^-3",
                                         species=:Ar, allowed_fail_count=1,
                                         initial_step_size=0.5, step_size=1.7)
        df, lso, p = @test_logs (:warn, r"Reset step size") match_mode=:any ideal_gas_referenced_nested_sampling(
            ls, p, 5, MCAtomGrandCanonicalMoves(), ac_save("stall"); initialize=false, stop_on_stall=false,
            record_move_rates=true)
        ac_clean("stall")
        @test nrow(df) == 5
        @test df.step_size[1] == 0.5    # the reset precedes the push
    end

    @testset "fail_count counts consecutive unmoved walks: $name" for (name, mk, mkp, routine) in cases
        # a kept clone's walk adds one and a tie eviction, which runs no walk, adds nothing;
        # five steps from six distinct keys cannot reach an all-tied live set
        Random.seed!(7108)
        ls = mk()
        p = mkp()
        expected = 0
        for _ in 1:5
            n = length(ls.walkers)
            out = step!(ls, p, routine)
            expected += (routine isa MCMixedMoves || length(ls.walkers) == n)   # a tie eviction changes the live count
            @test p.fail_count == expected
        end
    end

    @testset "a walk that moves resets fail_count" begin
        Random.seed!(7109)
        ls = ac_x_ls([1.5 + 1.1k for k in 1:6])
        p = nsp(50)
        p.fail_count = 5
        step!(ls, p, MCRandomWalkClone())
        @test p.fail_count == 0
        ls = ac_field_ls([[1.2 + 0.7j for j in 1:k] for k in 1:6])
        p = igp(50)
        p.fail_count = 5
        step!(ls, p, MCAtomGrandCanonicalMoves())
        @test p.fail_count == 0
    end

    @testset "the ideal-gas-referenced keep rule compares Ω = E - μN: μ = $mu" for mu in (0.1u"eV", -0.1u"eV")
        # with μ > 0 every Ω lies below every energy, with μ < 0 above it, so a rule that
        # read the energy, or the culled walker's energy as the ceiling, keeps nothing
        Random.seed!(7110)
        heights = [[1.2 + 0.7j for j in 1:k] for k in 1:6]
        omega(w) = w.energy - mu * w.list_num_par[1]
        ls = ac_field_ls(heights)
        sort!(ls.walkers, by=omega, rev=true)
        survivors = [snap(w) for w in ls.walkers[2:end]]
        p = igp(0; chemical_potential=mu)
        out = step!(ls, p, MCAtomGrandCanonicalMoves())
        @test out[1] isa Int && out[end] == -1 / 6
        @test any(s -> is_copy_of(ls.walkers[end], s), survivors)
        # the refill: the top two walkers are bit-exact duplicates
        top = mu > 0u"eV" ? 1 : 6
        rest = [k for k in 1:6 if k != top][1:4]
        ls = ac_field_ls([heights[top], heights[top], heights[rest]...])
        p = igp(0; chemical_potential=mu)
        step!(ls, p, MCAtomGrandCanonicalMoves())
        step!(ls, p, MCAtomGrandCanonicalMoves())
        @test length(ls.walkers) == 6
    end

    @testset "the charged prior mass of a one-move descent against its closed form (seeds 1 to 32)" begin
        devs = [ac_copy_dev(s) for s in 1:32]
        @test abs(mean(devs)) < 3.5
    end

    @testset "reference-fugacity descent in a field with short walks (seeds 1 to 32)" begin
        devs = [ac_igref_dev(s) for s in 1:32]
        @test abs(mean(first.(devs))) < 0.09
        @test abs(mean(last.(devs))) < 0.13
    end
end
