# One-body external fields on the atomistic grand-canonical and canonical paths
# (`ExternalFieldPotential`, `ZeroField`, `TabulatedPlanarField`, `TabulatedRadialField`).
# The evaluator unit tests and the field-first short-circuit live in test-EnergyEval.jl,
# which also defines the `NaNAboveField` fixture used below; the stream-independence
# fixture for the field lives in test-gc-stream-independence.jl.
#
# Calibration ledger (gates ship at >= 3x the maximum deviation over seeds {1, 2, 3}
# and the shipped seed, stated per statistic and per grid point):
# - Planar slit, reference-fugacity descent (K = 256, z0 V_acc = 4, 20-move walks, ideal
#   gas in U = 0.005 (z - 1) eV/Å on [1, 9] Å; the reduction at activities 0.7, 1 and 1.3
#   times z0 and T = 300 and 450 K, six grid points; seeds {1, 2, 3, 29001}): max |dev|
#   over the grid logXi 0.074, mean_N 0.148, var_N 0.224, mean_U 0.0022 eV; gates ship at
#   0.23, 0.45, 0.68 and 0.0065 eV.
# - Hard cylinder (K = 512, z0 V_acc = 4, all energies zero, an all-tied live set at the
#   first step; seeds {1, 2, 3, 29002}): max |dev| logXi 0.047, mean_N 0.258, var_N 0.843,
#   live-set mean occupancy against z0 V_acc 0.125; gates ship at 0.15, 0.78, 2.6 and 0.38.
#   Passing V_cell for V_acc would move logXi by z0 (V_cell - V_acc) = 3.3 at the reference
#   activity, so the gates bite on the convention they guard.
# - Radial wall (U = 0.002 r^2 eV tabulated at 1 Å nodes, r_end = 5 Å, K = 256, z0 V_acc
#   = 4; seeds {1, 2, 3, 29003}): max |dev| logXi 0.161, mean_N 0.198, var_N 0.348,
#   mean_U 0.0023 eV; gates ship at 0.49, 0.6, 1.05 and 0.0069 eV.
# - Chemical-potential ordering on the slit (mu = -0.02 eV; seeds {1, 2, 3, 29004}): max
#   |dev| logXi 0.057, mean_N 0.143, mean_U 0.0016 eV; gates ship at 0.18, 0.43 and
#   0.005 eV.
# - µVT driver at <N> = 3 (2e4 equilibration and 2e5 sampling steps at interval 10,
#   T = 300 K): slit (seeds {1, 2, 3, 29011}) max |dev| mean_N 0.048, var_N 0.088; hard
#   cylinder (seeds {1, 2, 3, 29012}) 0.037 and 0.061; gates ship at 0.15, 0.27, 0.12
#   and 0.19. The activity volume is z V_cell (the kernel inserts uniformly in the cell).
# - Kernel stationarity on the restricted measure (5e4 records of 10-step walks, ceiling out
#   of reach, ideal gas in the hard cylinder at z0 V_acc = 4; seeds {1, 2, 3, 29070}): max
#   |<N> - 4| 0.0151, max total variation to Poisson(4) 0.0087; gates ship at 0.046 and
#   0.027. A kernel ratio at z0 V_acc would give <N> = 2.18.
# - Every other testset is an exact contract and needs no calibration. The calibration
#   extracts the fixture block below verbatim and reruns it per seed.
@testset "external-field potentials on the atomistic samplers" begin
    using Random

    # >>> fixtures (the calibration script extracts this block verbatim)
    ef_kb = 8.617333262e-5
    ef_mass = 39.948u"u"
    ef_lam(T) = ustrip(u"Å", FreeBird.AnalysisTools._thermal_wavelength(ef_mass, T * u"K"))
    ef_mu_for(z, T) = (ef_kb * T * (log(z) + 3 * log(ef_lam(T)))) * u"eV"   # z in Å^-3
    ef_Ts = [300.0, 450.0]
    ef_scales = [0.7, 1.0, 1.3]

    # the slit: a linear planar field U = a (z - 1) on [1, 9] Å, +Inf outside, in a
    # 12 x 12 x 10 Å cell periodic in-plane
    ef_slit_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    ef_slit_at = FastSystem(atomic_system([:Ar => [6.0, 6.0, 5.0]u"Å"], ef_slit_box, (true, true, false)))
    ef_a = 0.005                     # eV/Å
    ef_W = 8.0                       # Å
    ef_A = 144.0                     # Å^2
    ef_slit = TabulatedPlanarField(3, [1.0, 9.0]u"Å", [0.0, ef_a * ef_W]u"eV")
    # the closed forms: I(T) = ∫ e^{-βU} dr and J(T) = ∫ U e^{-βU} dr over the slit
    ef_slit_I(T) = (β = 1 / (ef_kb * T); ef_A * (1 - exp(-β * ef_a * ef_W)) / (β * ef_a))
    ef_slit_J(T) = (β = 1 / (ef_kb * T);
                    ef_A * (1 - exp(-β * ef_a * ef_W) * (1 + β * ef_a * ef_W)) / (β^2 * ef_a))

    # the cylinder: a radial field about the z axis through (6, 6), r_end = 5 Å, in a
    # 12 x 12 x 10 Å cell periodic along the axis only
    ef_cyl_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    ef_cyl_at = FastSystem(atomic_system([:Ar => [6.0, 6.0, 5.0]u"Å"], ef_cyl_box, (false, false, true)))
    ef_rnodes = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    ef_hard = TabulatedRadialField(3, [6.0, 6.0, 0.0]u"Å", ef_rnodes * u"Å", zeros(6) * u"eV")
    ef_wall = TabulatedRadialField(3, [6.0, 6.0, 0.0]u"Å", ef_rnodes * u"Å", 0.002 .* ef_rnodes .^ 2 * u"eV")
    ef_L = 10.0
    # composite Simpson on each linear segment of the same table (exact up to 1e-12)
    function ef_radial_integrals(T; n=400)
        β = 1 / (ef_kb * T)
        U(r) = ustrip(u"eV", external_energy(ef_wall, SVector(6.0 + r, 6.0, 0.0)u"Å", ef_cyl_at))
        I = 0.0; J = 0.0
        for k in 1:5
            a, b = ef_rnodes[k], ef_rnodes[k+1]
            h = (b - a) / n
            for m in 0:n
                r = a + m * h
                wgt = (m == 0 || m == n) ? 1.0 : (isodd(m) ? 4.0 : 2.0)
                I += h / 3 * wgt * 2π * r * exp(-β * U(r))
                J += h / 3 * wgt * 2π * r * U(r) * exp(-β * U(r))
            end
        end
        return ef_L * I, ef_L * J
    end

    ef_mkempty(seed_at) = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                                     empty(position(seed_at, :)), empty(species(seed_at, :)),
                                     empty(mass(seed_at, :)))
    ef_save(tag) = SaveEveryN(df_filename="_extfield_$(tag).csv",
                              wk_filename="_extfield_$(tag).traj.extxyz",
                              ls_filename="_extfield_$(tag).ls.extxyz",
                              n_traj=10^7, n_snap=10^7, n_info=10^7)
    ef_clean(tag) = for f in ["_extfield_$(tag).csv", "_extfield_$(tag).traj.extxyz",
                              "_extfield_$(tag).ls.extxyz"]
        rm(f, force=true)
    end

    # one reference-fugacity descent of an ideal gas in a field, at z0 V_acc = z0V_acc
    function ef_descent(seed, pot, seed_at; K=256, z0V_acc=4.0, mc_steps=20, n_steps=8000,
                        mu=0.0u"eV", tag="d")
        Random.seed!(seed)
        ls = GenericAtomWalkers([AtomWalker{1}(ef_mkempty(seed_at)) for _ in 1:K], pot)
        Vacc = accessible_volume(pot.field, seed_at)
        z0 = z0V_acc / Vacc
        params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps, reference_activity=z0,
                                              species=:Ar, allowed_fail_count=10,
                                              chemical_potential=mu)
        df, lso, _ = ideal_gas_referenced_nested_sampling(ls, params, n_steps,
                                                          MCAtomGrandCanonicalMoves(), ef_save(tag))
        ef_clean(tag)
        live_e = [ustrip(u"eV", w.energy) for w in lso.walkers]
        live_n = [w.list_num_par[1] for w in lso.walkers]
        return df, live_e, live_n, z0, Vacc, lso
    end

    # the reduction at activities s z0 for s in ef_scales and T in ef_Ts, against the
    # closed forms ln Ξ = z I, <N> = Var N = z I, <U> = z J
    function ef_closure(df, live_e, live_n, z0, Vacc, IJ)
        mus = [ef_mu_for(s * ustrip(u"Å^-3", z0), T) for s in ef_scales, T in ef_Ts]
        dev = Dict(:logXi => zeros(3, 2), :mean_N => zeros(3, 2), :var_N => zeros(3, 2),
                   :mean_U => zeros(3, 2))
        for (j, T) in enumerate(ef_Ts)
            st = gc_thermodynamic_stats_ideal_ref(df, Vacc, ef_mass, z0, mus[:, j], [T]u"K";
                                                  live_emax=live_e, live_numbers=live_n)
            I, J = IJ(T)
            for (i, s) in enumerate(ef_scales)
                z = s * ustrip(u"Å^-3", z0)
                dev[:logXi][i, j] = st.logXi[i, 1] - z * I
                dev[:mean_N][i, j] = st.mean_N[i, 1] - z * I
                dev[:var_N][i, j] = st.var_N[i, 1] - z * I
                dev[:mean_U][i, j] = st.mean_U[i, 1] - z * J
            end
        end
        return dev
    end

    # the µVT driver on the same field at activity volume zV_cell (the kernel inserts
    # uniformly in the cell), from one atom at a finite energy
    function ef_muvt(seed, pot, seed_at, zV_cell, T)
        params = MuVTMCParameters([T], [zV_cell]; equilibrium_steps=20_000,
                                  sampling_steps=200_000, sampling_interval=10, random_seed=seed)
        w = AtomWalker{1}(deepcopy(seed_at))
        w.energy = interacting_energy(w.configuration, pot, w.list_num_par, w.frozen)
        return monte_carlo_sampling(MCAtomGrandCanonicalMoves(), w, pot, params)
    end
    # <<< fixtures

    @testset "ideal gas in a planar field closes against the closed forms (seed 29001)" begin
        pot = ExternalFieldPotential(IdealGasParameters(), ef_slit)
        df, live_e, live_n, z0, Vacc, lso = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ef_descent(29001, pot, ef_slit_at; tag="slit")
        end
        @test Vacc ≈ ef_A * ef_W * u"Å^3"
        @test nrow(df) > 256
        @test all(isfinite, df.emax) && all(isfinite, live_e)
        @test all(accessible(ef_slit, position(w.configuration, i), w.configuration)
                  for w in lso.walkers for i in 1:w.list_num_par[1])
        dev = ef_closure(df, live_e, live_n, z0, Vacc, T -> (ef_slit_I(T), ef_slit_J(T)))
        @test maximum(abs, dev[:logXi]) < 0.23
        @test maximum(abs, dev[:mean_N]) < 0.45
        @test maximum(abs, dev[:var_N]) < 0.68
        @test maximum(abs, dev[:mean_U]) < 0.0065
    end

    @testset "ideal gas in a hard cylinder: the restricted initializer and V_acc (seed 29002)" begin
        pot = ExternalFieldPotential(IdealGasParameters(), ef_hard)
        z0V_acc = 4.0
        df, live_e, live_n, z0, Vacc, lso = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ef_descent(29002, pot, ef_cyl_at; K=512, mc_steps=10, n_steps=50, z0V_acc=z0V_acc, tag="hard")
        end
        @test Vacc ≈ π * 25.0 * ef_L * u"Å^3"
        @test nrow(df) == 0                     # every energy is exactly zero: an all-tied live set
        @test all(iszero, live_e)
        @test all(accessible(ef_hard, position(w.configuration, i), w.configuration)
                  for w in lso.walkers for i in 1:w.list_num_par[1])
        # at the reference activity the reference-mass prefactor and the tail cancel exactly
        st = gc_thermodynamic_stats_ideal_ref(df, Vacc, ef_mass, z0,
                                              [ef_mu_for(s * ustrip(u"Å^-3", z0), 300.0) for s in ef_scales],
                                              [300.0]u"K"; live_emax=live_e, live_numbers=live_n)
        @test abs(st.logXi[2, 1] - z0V_acc) < 1e-12
        dev = ef_closure(df, live_e, live_n, z0, Vacc, T -> (ustrip(u"Å^3", Vacc), 0.0))
        @test maximum(abs, dev[:logXi]) < 0.15
        @test maximum(abs, dev[:mean_N]) < 0.78
        @test maximum(abs, dev[:var_N]) < 2.6
        @test all(st.mean_U .== 0.0)
        # the live set's occupancy is the restricted reference law, Poisson(z0 V_acc)
        @test abs(sum(live_n) / length(live_n) - z0V_acc) < 0.38
    end

    @testset "ideal gas against a tabulated radial wall (seed 29003)" begin
        pot = ExternalFieldPotential(IdealGasParameters(), ef_wall)
        df, live_e, live_n, z0, Vacc, lso = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ef_descent(29003, pot, ef_cyl_at; tag="wall")
        end
        @test all(isfinite, df.emax) && all(isfinite, live_e)
        dev = ef_closure(df, live_e, live_n, z0, Vacc, ef_radial_integrals)
        @test maximum(abs, dev[:logXi]) < 0.49
        @test maximum(abs, dev[:mean_N]) < 0.6
        @test maximum(abs, dev[:var_N]) < 1.05
        @test maximum(abs, dev[:mean_U]) < 0.0069
    end

    @testset "chemical-potential ordering with a field (seed 29004)" begin
        # at mu < 0 every Ω = E - µN is positive and the descent ends on the empty atom
        pot = ExternalFieldPotential(IdealGasParameters(), ef_slit)
        df, live_e, live_n, z0, Vacc, lso = @test_logs (:warn, r"stop_on_stall") match_mode=:any begin
            ef_descent(29004, pot, ef_slit_at; mu=-0.02u"eV", tag="omega")
        end
        @test issorted(df.omega, rev=true)
        @test all(isfinite, df.omega)
        dev = ef_closure(df, live_e, live_n, z0, Vacc, T -> (ef_slit_I(T), ef_slit_J(T)))
        @test maximum(abs, dev[:logXi]) < 0.18
        @test maximum(abs, dev[:mean_N]) < 0.43
        @test maximum(abs, dev[:mean_U]) < 0.005
    end

    @testset "µVT driver in a field: <N> and Var N against z I (seeds 29011, 29012)" begin
        T = 300.0
        V_cell = 1440.0
        # the slit: activity volume z V_cell with z I = 3
        z = 3.0 / ef_slit_I(T)
        out = ef_muvt(29011, ExternalFieldPotential(IdealGasParameters(), ef_slit), ef_slit_at, z * V_cell, T)
        @test abs(out.mean_N[1] - 3.0) < 0.15
        @test abs(out.var_N[1] - 3.0) < 0.27
        # the hard cylinder: z V_acc = 3
        z = 3.0 / (π * 25.0 * ef_L)
        out = ef_muvt(29012, ExternalFieldPotential(IdealGasParameters(), ef_hard), ef_cyl_at, z * V_cell, T)
        @test abs(out.mean_N[1] - 3.0) < 0.12
        @test abs(out.var_N[1] - 3.0) < 0.19
        @test all(out.mean_U .== 0.0)
    end

    @testset "zero-field weld: the shipped driver-pin fixtures, digit for digit" begin
        # the fixtures of test-atomistic-igref-driver-pins.jl, run once on the plain
        # Lennard-Jones potential and once wrapped with ZeroField, in the same process
        pin_box = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        pin_seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], pin_box, (true, true, true)))
        pin_V = 1728.0
        function pin_liveset(counts, pot)
            walkers = AtomWalker{1}[]
            for n in counts
                w = AtomWalker{1}(ef_mkempty(pin_seed_at))
                for _ in 1:n
                    pos = SVector(rand() * 12.0, rand() * 12.0, rand() * 12.0)u"Å"
                    FreeBird.AbstractWalkers.insert_particle!(w, pos, :Ar)
                end
                push!(walkers, w)
            end
            return GenericAtomWalkers(walkers, pot)
        end
        function pin_step_run(seed, counts, z0V, mc_steps, n_steps, wrap)
            Random.seed!(seed)
            lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
            ls = pin_liveset(counts, wrap ? ExternalFieldPotential(lj, ZeroField()) : lj)
            params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
                reference_activity=(z0V / pin_V)u"Å^-3", species=:Ar,
                allowed_fail_count=100_000, compression=:mean)
            out = (iters=Int[], emaxs=Float64[], npars=Int[], logts=Float64[])
            for k in 1:n_steps
                iter, emax, n_par, ls, params, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(
                    ls, params, MCAtomGrandCanonicalMoves(); ns_iteration=k, z0V=z0V)
                if !(iter isa Missing)
                    push!(out.iters, iter); push!(out.emaxs, ustrip(u"eV", emax))
                    push!(out.npars, n_par); push!(out.logts, log_t)
                end
            end
            return out, sort([ustrip(u"eV", w.energy) for w in ls.walkers])
        end
        for (seed, counts, z0V, mc, ns) in ((424261, [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8], 6.0, 60, 120),
                                             (424262, collect(1:12), 12.0, 40, 80))
            plain, live_p = pin_step_run(seed, counts, z0V, mc, ns, false)
            wrapped, live_w = pin_step_run(seed, counts, z0V, mc, ns, true)
            @test length(plain.iters) == ns
            @test wrapped.iters == plain.iters
            @test wrapped.npars == plain.npars
            @test wrapped.logts == plain.logts
            @test wrapped.emaxs == plain.emaxs
            @test live_w == live_p
        end
        # the kernel fixture: identical counters, acceptance and final count
        function pin_kernel_run(seed, wrap)
            Random.seed!(seed)
            w = AtomWalker{1}(ef_mkempty(pin_seed_at))
            for pos in ([2.0, 2.0, 2.0], [7.0, 7.0, 7.0], [2.0, 7.0, 2.0])
                FreeBird.AbstractWalkers.insert_particle!(w, SVector{3}(pos)u"Å", :Ar)
            end
            w.energy = 0.0u"eV"
            lj0 = LJParameters(epsilon=0.0)
            pot = wrap ? ExternalFieldPotential(lj0, ZeroField()) : lj0
            accept, rate, w, stats = MC_grand_canonical_walk!(500, w, pot, 1.0u"eV";
                z0V=4.0, species=:Ar, p_move=0.4, p_insert=0.3, step_size=0.8)
            return accept, rate, w.list_num_par[1], stats
        end
        @test pin_kernel_run(424253, true) == pin_kernel_run(424253, false)
        # the canonical walk: the wrapper's method reproduces the Lennard-Jones method
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
        function canon(pot)
            Random.seed!(424270)
            w = AtomWalker{1}(ef_mkempty(pin_seed_at))
            for i in 1:6
                FreeBird.AbstractWalkers.insert_particle!(w, SVector(2.0 * i, 3.0 + i, 6.0)u"Å", :Ar)
            end
            w.energy = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen)
            acc, rate, w = MC_random_walk!(400, w, pot, 0.6, w.energy + 0.02u"eV")
            return acc, rate, w.energy, position(w.configuration, :)
        end
        @test canon(ExternalFieldPotential(lj, ZeroField())) == canon(lj)
        # the grand-canonical initializer: a whole-cell field draws the base initializer's
        # stream (counts, positions, energies and the next draw), bounded and unbounded
        for (z0V, nmax) in ((6.0, typemax(Int64)), (12.0, 8))
            p = AtomisticIGRefGCNSParameters(reference_activity=(z0V / pin_V)u"Å^-3", species=:Ar, n_max=nmax)
            outs = map((lj, ExternalFieldPotential(lj, ZeroField()))) do pot
                Random.seed!(29050)
                ls = GenericAtomWalkers([AtomWalker{1}(ef_mkempty(pin_seed_at)) for _ in 1:8], pot)
                FreeBird.SamplingSchemes._init_atomistic_igref_walkers!(ls, p, z0V)
                ([w.list_num_par[1] for w in ls.walkers], [position(w.configuration, :) for w in ls.walkers],
                 [w.energy for w in ls.walkers], rand())
            end
            @test outs[1] == outs[2]
        end
    end

    @testset "canonical walk: incremental energies match a full recompute (seed 29020)" begin
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.0)
        pot = ExternalFieldPotential(lj, ef_wall)
        function mk()
            w = AtomWalker{1}(ef_mkempty(ef_cyl_at))
            for i in 1:6
                θ = 2π * i / 6
                FreeBird.AbstractWalkers.insert_particle!(w, SVector(6.0 + 3.0 * cos(θ), 6.0 + 3.0 * sin(θ), 1.5 * i)u"Å", :Ar)
            end
            w.energy = interacting_energy(w.configuration, pot, w.list_num_par, w.frozen)
            return w
        end
        w_inc = mk(); w_full = mk()
        emax = w_inc.energy + 0.05u"eV"
        Random.seed!(29020)
        acc_i, rate_i, w_inc = MC_random_walk!(2000, w_inc, pot, 0.8, emax)
        Random.seed!(29020)
        acc_f, rate_f, w_full = invoke(MC_random_walk!,
            Tuple{Int, AtomWalker{1}, FreeBird.AbstractPotentials.AbstractPotential, Float64, typeof(0.0u"eV")},
            2000, w_full, pot, 0.8, emax)
        @test acc_i == acc_f
        @test rate_i == rate_f
        @test 0.0 < rate_i < 1.0                  # both outcomes occur: moves out of the wall are rejected
        @test position(w_inc.configuration, :) == position(w_full.configuration, :)
        @test isapprox(ustrip(u"eV", w_inc.energy), ustrip(u"eV", w_full.energy); rtol=1e-12)
        @test all(accessible(ef_wall, position(w_inc.configuration, i), w_inc.configuration) for i in 1:6)
    end

    @testset "no infinities: hard regions and NaN-returning fields" begin
        # a seeded hard-cylinder descent of an interacting gas: no ledger row and no live
        # walker carries ±Inf or NaN, and every atom stays in the region
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.0)
        pot = ExternalFieldPotential(lj, ef_hard)
        # the save layer's zero-atom warning fires or not by trajectory (so by Julia version):
        # swallowed without a pattern
        df, live_e, live_n, z0, Vacc, lso = @test_logs match_mode=:any begin
            ef_descent(29030, pot, ef_cyl_at; K=64, z0V_acc=6.0, mc_steps=40, n_steps=400, tag="inf")
        end
        @test nrow(df) > 100
        @test all(isfinite, df.emax)
        @test all(isfinite, live_e)
        @test all(accessible(ef_hard, position(w.configuration, i), w.configuration)
                  for w in lso.walkers for i in 1:w.list_num_par[1])
        # kernel level: a region so small that every insertion lands outside it and every
        # displacement leaves it; the walker's energy is unchanged bit for bit
        tiny = TabulatedRadialField(3, [6.0, 6.0, 0.0]u"Å", [0.0, 0.01]u"Å", [0.0, 0.0]u"eV")
        tpot = ExternalFieldPotential(IdealGasParameters(), tiny)
        Random.seed!(29031)
        w = AtomWalker{1}(ef_mkempty(ef_cyl_at))
        FreeBird.AbstractWalkers.insert_particle!(w, SVector(6.0, 6.0, 5.0)u"Å", :Ar)
        w.energy = 0.0u"eV"
        e0 = w.energy
        # z0V = 1e6 puts every insertion ratio above one (a rejection is the ceiling's) and
        # every deletion ratio near 5e-5; an open deletion channel keeps the ratios nonzero
        _, _, w, stats = MC_grand_canonical_walk!(300, w, tpot, 1.0u"eV"; z0V=1.0e6, species=:Ar,
                                                  p_move=0.5, p_insert=0.49, step_size=0.5)
        @test stats.delete_accepted == 0
        @test stats.insert_attempted > 0 && stats.insert_accepted == 0
        @test stats.move_attempted > 0 && stats.move_accepted == 0
        @test w.energy === e0
        @test position(w.configuration, 1) == SVector(6.0, 6.0, 5.0)u"Å"
        # a NaN-returning field: the ceiling tests reject NaN proposals (the grand-canonical
        # kernel, the wrapper's canonical walk and the generic full-recompute walk)
        npot = ExternalFieldPotential(IdealGasParameters(), NaNAboveField(6.0u"Å"))
        box12 = [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]]u"Å"
        nan_at = FastSystem(atomic_system([:Ar => [3.0, 6.0, 6.0]u"Å", :Ar => [4.0, 2.0, 9.0]u"Å"], box12, (true, true, true)))
        Random.seed!(29032)
        w = AtomWalker{1}(deepcopy(nan_at)); w.energy = 0.0u"eV"
        _, _, w, stats = MC_grand_canonical_walk!(400, w, npot, 1.0u"eV"; z0V=1.0e6, species=:Ar,
                                                  p_move=0.5, p_insert=0.49, step_size=1.0)
        @test stats.insert_accepted > 0 && stats.insert_accepted < stats.insert_attempted
        @test w.energy == 0.0u"eV"
        @test all(position(w.configuration, i)[1] <= 6.0u"Å" for i in 1:length(w.configuration))
        for walk in (MC_random_walk!,
                     (n, w, p, s, e) -> invoke(MC_random_walk!,
                         Tuple{Int, AtomWalker{1}, FreeBird.AbstractPotentials.AbstractPotential, Float64, typeof(0.0u"eV")},
                         n, w, p, s, e))
            Random.seed!(29033)
            w = AtomWalker{1}(deepcopy(nan_at)); w.energy = 0.0u"eV"
            acc, rate, w = walk(400, w, npot, 1.0, 1.0u"eV")
            @test 0.0 < rate < 1.0
            @test w.energy == 0.0u"eV"
            @test all(position(w.configuration, i)[1] <= 6.0u"Å" for i in 1:2)
        end
    end

    @testset "routine and warning contracts" begin
        # the Galilean burst is refused with a field, at the driver and at the step
        pot = ExternalFieldPotential(IdealGasParameters(), ef_hard)
        ls = GenericAtomWalkers([AtomWalker{1}(ef_mkempty(ef_cyl_at)) for _ in 1:4], pot)
        params = AtomisticIGRefGCNSParameters(mc_steps=5, reference_activity=0.004u"Å^-3", species=:Ar)
        gal = MCAtomGrandCanonicalMoves(galilean_steps=3)
        @test_throws ArgumentError ideal_gas_referenced_nested_sampling(ls, params, 2, gal, ef_save("gal"))
        @test_throws ArgumentError FreeBird.SamplingSchemes.nested_sampling_step!(ls, params, gal)
        ef_clean("gal")
        # the bounded construction's guard reads the restricted law (z0 V_acc), not the cell's:
        # a narrow cylinder with z0 V_acc = 2 and n_max = 5 is workable (P(N <= 5) = 0.983),
        # though z0 V_cell = 92 would put the whole cell's law far above the cap
        narrow = TabulatedRadialField(3, [6.0, 6.0, 0.0]u"Å", [0.0, 1.0]u"Å", [0.0, 0.0]u"eV")
        Vn = ustrip(u"Å^3", accessible_volume(narrow, ef_cyl_at))
        npot = ExternalFieldPotential(IdealGasParameters(), narrow)
        Random.seed!(29060)
        ls = GenericAtomWalkers([AtomWalker{1}(ef_mkempty(ef_cyl_at)) for _ in 1:8], npot)
        pb = AtomisticIGRefGCNSParameters(mc_steps=5, reference_activity=(2.0 / Vn)u"Å^-3",
                                          species=:Ar, n_max=5, allowed_fail_count=3)
        df, lso, _ = @test_logs match_mode=:any ideal_gas_referenced_nested_sampling(
            ls, pb, 5, MCAtomGrandCanonicalMoves(), ef_save("nmax"))
        ef_clean("nmax")
        @test all(w.list_num_par[1] <= 5 for w in lso.walkers)
        pz = AtomisticIGRefGCNSParameters(mc_steps=5, reference_activity=(25.0 / Vn)u"Å^-3",
                                          species=:Ar, n_max=1)       # P(N <= 1) = 3.6e-10
        err = try
            ideal_gas_referenced_nested_sampling(ls, pz, 2, MCAtomGrandCanonicalMoves(), ef_save("nmax0"))
            nothing
        catch e
            e
        end
        ef_clean("nmax0")
        @test err isa ArgumentError && occursin("z0V_acc", err.msg)
        # the minimum-image check reads periodic axes only
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)              # range 6.25 Å
        wide = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"],
            [[16.0, 0.0, 0.0], [0.0, 16.0, 0.0], [0.0, 0.0, 10.0]]u"Å", (true, true, false)))
        @test_logs FreeBird.SamplingSchemes._warn_min_image_cutoff(lj, wide)
        closed = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"],
            [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]]u"Å", (false, false, false)))
        @test_logs FreeBird.SamplingSchemes._warn_min_image_cutoff(lj, closed)
        narrow = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"],
            [[16.0, 0.0, 0.0], [0.0, 16.0, 0.0], [0.0, 0.0, 10.0]]u"Å", (true, true, true)))
        @test_logs (:warn, r"minimum-image") FreeBird.SamplingSchemes._warn_min_image_cutoff(lj, narrow)
        @test_logs (:warn, r"minimum-image") FreeBird.SamplingSchemes._warn_min_image_cutoff(
            ExternalFieldPotential(lj, ZeroField()), narrow)
    end

    @testset "absolute pins of the field path (captured on this change set)" begin
        # >>> pins (the capture script extracts this block verbatim)
        pin_lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.0)
        function field_liveset(counts, pot, seed_at, draw)
            walkers = AtomWalker{1}[]
            for n in counts
                w = AtomWalker{1}(ef_mkempty(seed_at))
                for _ in 1:n
                    FreeBird.AbstractWalkers.insert_particle!(w, draw(), :Ar)
                end
                push!(walkers, w)
            end
            return GenericAtomWalkers(walkers, pot)
        end
        slit_draw() = SVector(rand() * 12.0, rand() * 12.0, 1.0 + rand() * 8.0)u"Å"
        cyl_draw() = (ρ = 5.0 * sqrt(rand()); θ = 2π * rand();
                      SVector(6.0 + ρ * cos(θ), 6.0 + ρ * sin(θ), rand() * 10.0)u"Å")
        function field_step_run(seed, pot, seed_at, draw, counts, z0V, mc_steps, n_steps)
            Random.seed!(seed)
            ls = field_liveset(counts, pot, seed_at, draw)
            V_cell = ustrip(u"Å^3", accessible_volume(ZeroField(), seed_at))
            params = AtomisticIGRefGCNSParameters(mc_steps=mc_steps,
                reference_activity=(z0V / V_cell)u"Å^-3", species=:Ar, allowed_fail_count=100_000)
            iters = Int[]; emaxs = Float64[]; npars = Int[]
            for k in 1:n_steps
                iter, emax, n_par, ls, params, log_t = FreeBird.SamplingSchemes.nested_sampling_step!(
                    ls, params, MCAtomGrandCanonicalMoves(); ns_iteration=k, z0V=z0V)
                if !(iter isa Missing)
                    push!(iters, iter); push!(emaxs, ustrip(u"eV", emax)); push!(npars, n_par)
                end
            end
            return iters, emaxs, npars, sort([ustrip(u"eV", w.energy) for w in ls.walkers])
        end
        function field_kernel_run(seed, pot, seed_at)
            Random.seed!(seed)
            w = AtomWalker{1}(ef_mkempty(seed_at))
            for pos in ([6.0, 6.0, 2.0], [8.0, 6.0, 5.0], [6.0, 4.0, 7.0])
                FreeBird.AbstractWalkers.insert_particle!(w, SVector{3}(pos)u"Å", :Ar)
            end
            w.energy = interacting_energy(w.configuration, pot, w.list_num_par, w.frozen)
            accept, rate, w, stats = MC_grand_canonical_walk!(500, w, pot, w.energy + 0.05u"eV";
                z0V=6.0, species=:Ar, p_move=0.4, p_insert=0.3, step_size=0.8)
            return accept, rate, w.list_num_par[1], stats
        end
        pin_slit() = field_step_run(29041, ExternalFieldPotential(pin_lj, ef_slit), ef_slit_at, slit_draw,
                                    collect(1:12), 8.0, 40, 60)
        pin_wall() = field_step_run(29042, ExternalFieldPotential(pin_lj, ef_wall), ef_cyl_at, cyl_draw,
                                    collect(1:10), 6.0, 40, 60)
        pin_kern() = field_kernel_run(29043, ExternalFieldPotential(pin_lj, ef_wall), ef_cyl_at)
        # <<< pins
        PIN_SLIT_N = 60
        PIN_SLIT_NPAR = [
            9, 11, 12, 4, 8, 7, 10, 13, 7, 6, 6, 8, 10, 8, 7, 9, 4, 5, 6, 3, 7, 5, 5, 5, 5, 5, 3, 3, 3, 4,
            2, 3, 3, 3, 2, 2, 3, 2, 2, 3, 2, 3, 3, 2, 1, 2, 2, 3, 3, 4, 1, 2, 1, 1, 2, 2, 2, 1, 3, 3]
        PIN_SLIT_EMAX = [
            146.05173074288888, 73.6745214906074, 48.34824542369292, 5.712463847140378, 2.795579780143612,
            0.6420517802068397, 0.62913953881569, 0.5532814665460773, 0.44512739853857847,
            0.3460473284805831, 0.1847404815013459, 0.15327264746598315, 0.1514945438309219,
            0.13518064573942612, 0.12789542346072533, 0.1251225457724216, 0.11590804455492748,
            0.11440787458400413, 0.1139538155008813, 0.10799997334552405, 0.10466156121262554,
            0.09848883544253989, 0.08740039050988971, 0.08676468644505324, 0.08467282794961686,
            0.07836979592572417, 0.07092602108792201, 0.057429582908069154, 0.054538341584756296,
            0.05100159607010206, 0.04905043812997657, 0.04179695734065998, 0.041598149228332805,
            0.03881965977568919, 0.03685719158895083, 0.03582598897103301, 0.03343896277776327,
            0.03310127088528072, 0.032854789850332286, 0.03227613459430116, 0.032196434601111774,
            0.02821914414039305, 0.028197202640980842, 0.027138166429107147, 0.025900552098307536,
            0.02500747857959357, 0.024945829926076565, 0.022524450874416367, 0.021371516598961953,
            0.021102172600227898, 0.0190228180291609, 0.018713900592657654, 0.018584400193994703,
            0.017403488511678313, 0.016781980008038662, 0.01555607987087974, 0.014499850057240055,
            0.010993913998919025, 0.010781314121557023, 0.010492366578840655]
        PIN_SLIT_LIVE = [
            -0.0039874755915961435, 0.0, 0.0004626295943304798, 0.0004763513691387944,
            0.0005638737248862737, 0.0006035206310032671, 0.0020770513378906667, 0.005143714494497007,
            0.006123024682195681, 0.008687850073374989, 0.009170037274932702, 0.009297625402276168]
        PIN_WALL_N = 43
        PIN_WALL_NPAR = [
            5, 9, 8, 7, 6, 4, 10, 7, 7, 4, 5, 7, 2, 4, 3, 2, 3, 1, 2, 2, 1, 2, 5, 1, 3, 4, 1, 2, 1, 1, 1, 1,
            1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1]
        PIN_WALL_EMAX = [
            634967.2787520837, 591.1752049400078, 482.0540661958828, 109.09965332334616, 42.81149327804999,
            19.974476989558557, 0.3913391336736295, 0.2874234651645596, 0.18149149381848123,
            0.12459611762241478, 0.11757838488931746, 0.11011039150865856, 0.07424959110448404,
            0.07018992734701728, 0.06289211740239992, 0.05466039578094743, 0.04985789593150186,
            0.04937299602332831, 0.04804839891581934, 0.04456522603044666, 0.03613485222111371,
            0.033816300882912737, 0.03319924197774608, 0.032137111641218574, 0.027199616841648688,
            0.02567251482297891, 0.0250599523713161, 0.024996076814299408, 0.02347299238939422,
            0.022801255126998557, 0.022411736463821164, 0.02005989762639165, 0.016869637785457376,
            0.014637545180953298, 0.014331614935435362, 0.014172475253665742, 0.012133448184368614,
            0.009554179324218353, 0.008752538682666289, 0.007279564524890998, 0.003935063955954682,
            0.003267575308172903, 0.0015817993581199936]
        PIN_WALL_LIVE = [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        PIN_KERN = (true, 0.522, 2)
        PIN_KERN_STATS = (move_attempted = 204, move_accepted = 174, insert_attempted = 155, insert_accepted = 43, insert_biased_attempted = 0, insert_biased_accepted = 0, delete_attempted = 131, delete_accepted = 44)

        # every Julia version: same-process replay identity and the version-robust contracts
        for f in (pin_slit, pin_wall)
            a = f(); b = f()
            @test a == b
            iters, emaxs, npars, live = a
            @test iters == collect(1:length(iters))
            @test issorted(emaxs, rev=true)
            @test all(isfinite, emaxs) && all(isfinite, live)
        end
        @test pin_kern() == pin_kern()
        # captured on Julia 1.10 (reproduced digit-identically in two separate processes);
        # accept/reject cascades differ on later versions, so the trajectory vectors are
        # guarded (the policy of the shipped pin files): smooth accumulations at rtol 1e-12,
        # integer sequences and counters exact
        if VERSION < v"1.11"
            for (f, N, npar, emax, live) in ((pin_slit, PIN_SLIT_N, PIN_SLIT_NPAR, PIN_SLIT_EMAX, PIN_SLIT_LIVE),
                                             (pin_wall, PIN_WALL_N, PIN_WALL_NPAR, PIN_WALL_EMAX, PIN_WALL_LIVE))
                iters, emaxs, npars, lv = f()
                @test length(iters) == N
                @test npars == npar
                @test all(isapprox.(emaxs, emax; rtol=1e-12, atol=0.0))
                @test all(isapprox.(lv, live; rtol=1e-12, atol=0.0))
            end
            accept, rate, n, stats = pin_kern()
            @test (accept, rate, n) == PIN_KERN
            for (k, v) in pairs(PIN_KERN_STATS)
                @test getproperty(stats, k) == v
            end
        end
    end

    @testset "kernel stationarity on the restricted measure (seed 29070)" begin
        # with the ceiling out of reach, the grand-canonical kernel's stationary law on an ideal
        # gas in the hard cylinder is the restricted reference law Poisson(z0 V_acc) when its
        # insertion ratio uses z0 V_cell; a ratio at z0 V_acc would give <N> = z0 V_acc^2 / V_cell
        pot = ExternalFieldPotential(IdealGasParameters(), ef_hard)
        Vacc = ustrip(u"Å^3", accessible_volume(ef_hard, ef_cyl_at))
        Vcell = 1440.0
        z0V_acc = 4.0
        Random.seed!(29070)
        w = AtomWalker{1}(deepcopy(ef_cyl_at))
        w.energy = 0.0u"eV"
        n_rec = 50_000
        counts = zeros(Int, 30)
        for _ in 1:n_rec
            MC_grand_canonical_walk!(10, w, pot, 1.0e6u"eV"; z0V=z0V_acc * Vcell / Vacc, species=:Ar,
                                     p_move=0.4, p_insert=0.3, step_size=1.0)
            counts[min(w.list_num_par[1], 29) + 1] += 1
        end
        p = counts ./ n_rec
        pois = [exp(n * log(z0V_acc) - z0V_acc - sum(log, 1:n; init=0.0)) for n in 0:29]
        @test abs(sum((0:29) .* p) - z0V_acc) < 0.046
        @test 0.5 * sum(abs.(p .- pois)) < 0.027
        @test all(accessible(ef_hard, position(w.configuration, i), w.configuration) for i in 1:w.list_num_par[1])
    end
end
