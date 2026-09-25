# Fixed-temperature grand-canonical (muVT) Metropolis sampling driver.
#
# Calibration ledger (gates at >= 3x the max three-seed deviation, per stat):
# - Ideal-gas closure, T = [300, 400] K, zV = [3, 5], seeds 95001/95002/95003:
#   max devs mean 0.0411, var 0.1479; gates 0.13 and 0.45.
# - Dilute-LJ second-virial mean gate (epsilon = 0.01 eV, sigma = 2.0 A,
#   cutoff 2.5 so r_c = 5.0 A = L/2, guard-clean; b2 = +1.2458 A^3 at 300 K,
#   N_ref = 0.80159 at zV = 0.8), seeds 96001/96002/96003: max dev 0.00968;
#   gate 0.03.
struct LatticeMuVTFieldHamiltonian <: ClassicalHamiltonian
    epsilon::Float64
end

FreeBird.EnergyEval.interacting_energy(
    lattice::AtomicLattice, h::LatticeMuVTFieldHamiltonian) =
    h.epsilon * n_occupied(lattice) * u"eV"

@testset "muVT Metropolis sampling driver" begin
    using Random

    box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    pbc = (true, true, true)
    mk1() = AtomWalker(FastSystem(atomic_system([:Ar => [5.0, 5.0, 5.0]u"Å"], box, pbc)))
    lj0 = LJParameters(epsilon=0.0)

    @testset "parameter and entry validation" begin
        @test_throws ArgumentError MuVTMCParameters([300.0], [1.0, 2.0])
        @test_throws ArgumentError MuVTMCParameters([300.0], [-1.0])
        @test_throws ArgumentError MuVTMCParameters([300.0], [1.0]; sampling_interval=0)
        params = MuVTMCParameters([300.0], [1.0])
        # omitting the routine hits the muVT catch-all, not a bare MethodError
        @test_throws ArgumentError monte_carlo_sampling(mk1(), lj0, params)
        # an empty starting walker carries no species record
        empty_sys = FastSystem(cell_vectors(mk1().configuration), pbc,
                               empty(position(mk1().configuration, :)),
                               empty(species(mk1().configuration, :)),
                               empty(mass(mk1().configuration, :)))
        w_empty = AtomWalker{1}(empty_sys)
        @test_throws ArgumentError monte_carlo_sampling(MCAtomGrandCanonicalMoves(),
            w_empty, lj0, params)
    end

    @testset "ideal-gas closure over a temperature ladder (seeded, calibrated)" begin
        params = MuVTMCParameters([300.0, 400.0], [3.0, 5.0];
            equilibrium_steps=20_000, sampling_steps=200_000, sampling_interval=10,
            random_seed=95001)
        w = mk1(); w.energy = 0.0u"eV"
        out = monte_carlo_sampling(MCAtomGrandCanonicalMoves(), w, lj0, params)
        for (i, zV) in enumerate([3.0, 5.0])
            @test abs(out.mean_N[i] - zV) < 0.13
            @test abs(out.var_N[i] - zV) < 0.45
            @test abs(sum(out.p_N[i]) - 1.0) < 1e-12
            @test length(out.N_series[i]) == 200_000 ÷ 10
            @test out.N_support[i] == 0:maximum(out.N_series[i])
        end
        @test all(out.mean_U .== 0.0)
    end

    @testset "dilute-LJ second-virial mean gate (seeded, calibrated)" begin
        lj = LJParameters(epsilon=0.01, sigma=2.0, cutoff=2.5, shift=true)
        kb = 8.617333262e-5
        T = 300.0
        β = 1 / (kb * T)
        u_of(r) = ustrip(u"eV", FreeBird.AbstractPotentials.lj_energy((r)u"Å", lj))
        simpson(f, a, b; n=4000) = (h = (b - a) / n;
            h / 3 * sum(k -> f(a + k * h) * (k == 0 || k == n ? 1 : (isodd(k) ? 4 : 2)), 0:n))
        b2 = 2π * simpson(r -> (exp(-β * u_of(r)) - 1) * r^2, 1.0e-6, 5.0)
        zV = 0.8
        N_ref = zV + 2 * zV^2 * b2 / 1000.0
        params = MuVTMCParameters([T], [zV];
            equilibrium_steps=20_000, sampling_steps=400_000, sampling_interval=10,
            random_seed=96001)
        w = mk1(); w.energy = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen)
        out = monte_carlo_sampling(MCAtomGrandCanonicalMoves(), w, lj, params)
        @test abs(out.mean_N[1] - N_ref) < 0.03
        # production acceptance counters are populated and channel-consistent
        acc = out.acceptance[1]
        @test acc.insert_attempted > 0 && acc.delete_attempted > 0
        @test acc.insert_accepted <= acc.insert_attempted
        @test acc.delete_accepted <= acc.delete_attempted
        # the returned walker's incremental energy stays anchored
        wf = out.walkers[1]
        E_re = interacting_energy(wf.configuration, lj, wf.list_num_par, wf.frozen) +
               wf.energy_frozen_part
        @test abs(ustrip(u"eV", wf.energy - E_re)) < 5e-10
    end

    @testset "phase seeding and frozen production (contract emulation)" begin
        # The driver's documented contract: equilibration seeded random_seed in
        # ten adaptation blocks on the displacement-only rate; production seeded
        # random_seed + 1 with the step size frozen. An independent emulation of
        # that contract must reproduce the recorded series exactly.
        lj = LJParameters(epsilon=0.01, sigma=2.0, cutoff=2.5, shift=true)
        seed = 97001
        params = MuVTMCParameters([300.0], [2.0];
            equilibrium_steps=2_000, sampling_steps=8_000, sampling_interval=20,
            random_seed=seed)
        w = mk1(); w.energy = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen)
        out = monte_carlo_sampling(MCAtomGrandCanonicalMoves(), w, lj, params)

        emu_params = MuVTMCParameters([300.0], [2.0];
            equilibrium_steps=2_000, sampling_steps=8_000, sampling_interval=20,
            random_seed=seed)
        we = mk1(); we.energy = interacting_energy(we.configuration, lj, we.list_num_par, we.frozen)
        sp_e = species(we.configuration, 1)
        Random.seed!(seed)
        for _ in 1:10
            _, _, stats = MC_muVT_walk!(200, we, lj, 300.0; zV=2.0, species=sp_e,
                                        step_size=emu_params.step_size)
            if stats.move_attempted > 0
                FreeBird.SamplingSchemes.adjust_step_size(emu_params,
                    stats.move_accepted / stats.move_attempted; range=emu_params.accept_range)
            end
        end
        Random.seed!(seed + 1)
        ns_emu = Int[]
        for _ in 1:(8_000 ÷ 20)
            MC_muVT_walk!(20, we, lj, 300.0; zV=2.0, species=sp_e,
                          step_size=emu_params.step_size)
            push!(ns_emu, we.list_num_par[1])
        end
        @test out.N_series[1] == ns_emu
        @test out.mean_N[1] == mean(ns_emu)
    end

    @testset "same-seed determinism" begin
        lj = LJParameters(epsilon=0.01, sigma=2.0, cutoff=2.5, shift=true)
        outs = map(1:2) do _
            params = MuVTMCParameters([300.0, 350.0], [1.0, 2.0];
                equilibrium_steps=1_000, sampling_steps=5_000, sampling_interval=10,
                random_seed=98001)
            w = mk1(); w.energy = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen)
            monte_carlo_sampling(MCAtomGrandCanonicalMoves(), w, lj, params)
        end
        @test outs[1].mean_N == outs[2].mean_N
        @test outs[1].N_series == outs[2].N_series
        @test outs[1].acceptance == outs[2].acceptance
    end

    @testset "cutoff guard at driver entry" begin
        over = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5, shift=true)   # 6.25 > 5.0
        params = MuVTMCParameters([300.0], [1.0];
            equilibrium_steps=100, sampling_steps=200, sampling_interval=10)
        w = mk1(); w.energy = interacting_energy(w.configuration, over, w.list_num_par, w.frozen)
        @test_logs (:warn, r"minimum-image") match_mode = :any begin
            monte_carlo_sampling(MCAtomGrandCanonicalMoves(), w, over, params)
        end
    end
end


@testset "AtomicLattice muVT Metropolis sampling driver" begin
    lattice = AtomicLattice{1,SquareLattice}(
        lattice_atom="Pd",
        type_of_sites=["hollow"],
        coverage=0.25,
        adsorbate_atoms=["O"],
        supercell_dimensions=(2, 2, 1),
        lattice_constant=3.947,
        periodicity=(true, true, false),
        num_nearest_neighbors=2,
    )
    h = LatticeMuVTFieldHamiltonian(0.02)
    routine = MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25)

    @testset "parameter and option validation" begin
        params = MetropolisMCParameters(
            [600.0]; equilibrium_steps=100, sampling_steps=200,
            chemical_potentials=[-0.03, 0.0])
        @test params.chemical_potentials == [-0.03, 0.0]
        @test_throws ArgumentError MetropolisMCParameters(
            [600.0]; chemical_potentials=Float64[])
        @test_throws ArgumentError MetropolisMCParameters(
            [600.0]; chemical_potentials=[NaN])
        @test_throws ArgumentError monte_carlo_sampling(
            routine, lattice, h, MetropolisMCParameters([600.0]))
        @test_throws ArgumentError μvt_monte_carlo(
            MCGrandCanonicalMoves(p_bias=0.5), lattice, h, 600.0, 10, 1)
        @test_throws ArgumentError μvt_monte_carlo(
            MCGrandCanonicalMoves(clusters_freq=1), lattice, h, 600.0, 10, 1)
    end

    @testset "multi-species guard" begin
        multi = AtomicLattice{2,SquareLattice}(
            lattice_atom="Pd",
            type_of_sites=["hollow"],
            components=[1, 1],
            adsorbate_atoms=["O", "H"],
            supercell_dimensions=(2, 2, 1),
            lattice_constant=3.947,
            periodicity=(true, true, false),
            num_nearest_neighbors=2,
        )
        err = try
            μvt_monte_carlo(routine, multi, h, 600.0, 10, 1)
            nothing
        catch exception
            exception
        end
        @test err isa ArgumentError
        @test occursin("one-species AtomicLattice", sprint(showerror, err))

        params = MetropolisMCParameters(
            [600.0]; equilibrium_steps=0, sampling_steps=10,
            chemical_potentials=[-0.03])
        @test_throws ArgumentError monte_carlo_sampling(
            routine, multi, h, params)
        @test_throws ArgumentError FreeBird.SamplingSchemes.metropolis_hastings(
            routine, multi, 0.0, deepcopy(multi), h, 1.0, -0.03, :move)
    end

    @testset "reverse-proposal detailed balance" begin
        beta = 1 / (8.617_333_262e-5 * 600.0)
        p_insert = 0.2
        p_delete = 0.3
        n = 1
        n_sites = 4
        delta_omega = 0.017
        a_forward = FreeBird.SamplingSchemes._lattice_muvt_acceptance_probability(
            delta_omega, beta, :insert, n, n_sites, p_insert, p_delete)
        a_reverse = FreeBird.SamplingSchemes._lattice_muvt_acceptance_probability(
            -delta_omega, beta, :delete, n + 1, n_sites,
            p_insert, p_delete)
        q_forward = p_insert / (n_sites - n)
        q_reverse = p_delete / (n + 1)
        @test isapprox(q_forward * a_forward,
                       exp(-beta * delta_omega) * q_reverse * a_reverse;
                       rtol=1e-14)
    end

    @testset "thinned trajectory stores bare energy" begin
        energies, configs, coverages, accepted = μvt_monte_carlo(
            routine, lattice, h, 600.0, 23, 20260914;
            μ=-0.03, record_interval=5)
        @test length(energies) == length(configs) == length(coverages) == 5
        @test 0 <= accepted <= 23
        @test all(energies .== 0.02 .* n_occupied.(configs))
        @test all(coverages .== n_occupied.(configs) ./ num_sites(lattice))
    end

    @testset "Python calculator energy path" begin
        ase_lj = FreeBird.EnergyEval.pyimport(
            "ase.calculators.lj").LennardJones()
        calc = PyMLPotential(
            FreeBird.AbstractPotentials.ASEcalculator(ase_lj))
        energies, configs, coverages, accepted = μvt_monte_carlo(
            routine, lattice, calc, 600.0, 4, 20260916;
            μ=-0.03, record_interval=2)
        @test length(energies) == length(configs) == length(coverages) == 2
        @test all(isfinite, energies)
        @test 0 <= accepted <= 4
        @test all(!FreeBird.AbstractWalkers.pyis(
            lattice.ase_lattice, config.ase_lattice) for config in configs)
    end

    @testset "independent-site exact mean and high-level schema" begin
        temperature = 600.0
        mu = -0.03
        params = MetropolisMCParameters(
            [temperature]; equilibrium_steps=5_000, sampling_steps=50_000,
            chemical_potentials=[mu], random_seed=20260915)
        df, finals = monte_carlo_sampling(
            routine, lattice, h, params; sampling_interval=5)
        exact_coverage = 1 / (1 + exp((h.epsilon - mu) /
                                     (8.617_333_262e-5 * temperature)))
        @test nrow(df) == 1
        @test names(df) == ["temperature", "energy", "c_omega", "cov",
                            "acceptance_rate", "chemical_potential"]
        @test abs(df.cov[1] - exact_coverage) < 0.03
        @test isapprox(df.energy[1], h.epsilon * num_sites(lattice) * df.cov[1];
                       atol=1e-12)
        @test isfinite(df.c_omega[1]) && df.c_omega[1] >= 0.0
        @test 0.0 <= df.acceptance_rate[1] <= 1.0
        @test haskey(finals, mu)
    end
end
