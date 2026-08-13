@testset "nvt_monte_carlo.jl tests" begin

    @testset "MCMetropolisMCParameters struct tests" begin
        # Basic tests
        params = MetropolisMCParameters(
            [1.0, 2.0, 3.0],
            equilibrium_steps=100,
            sampling_steps=200,
            step_size=0.1,
            step_size_up=0.2,
            accept_range=(0.3, 0.4)
        )

        @test params.temperatures == [1.0, 2.0, 3.0]
        @test params.equilibrium_steps == 100
        @test params.sampling_steps == 200
        @test params.step_size == 0.1
        @test params.step_size_up == 0.2
        @test params.accept_range == (0.3, 0.4)
    end
    @testset "nvt_monte_carlo function tests" begin
        
        # Setup
        lattice = SLattice{SquareLattice}(
                supercell_dimensions=(2, 1, 1),
                components=[[true, false]]
            )
        
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

        mc_routine = MCNewSample()

        @testset "Basic functionality" begin

            energies, configs, accepted = nvt_monte_carlo(mc_routine, lattice, ham, 300.0, 100, 42)

            @test length(energies) == 100
            @test length(configs) == 100
            @test 0 ≤ accepted ≤ 100
            @test all(isfinite, energies)
        end
    
        @testset "Temperature effects" begin

            energies_cold, _, accepted_cold = nvt_monte_carlo(mc_routine, lattice, ham, 1.0, 50, 42)
            energies_hot, _, accepted_hot = nvt_monte_carlo(mc_routine, lattice, ham, 1000.0, 50, 42)

            @test std(energies_cold) ≤ std(energies_hot)
            @test accepted_cold ≤ accepted_hot
        end
    
        @testset "Conservation & Reproducibility" begin

            energies1, configs1, accepted1 = nvt_monte_carlo(mc_routine, lattice, ham, 300.0, 50, 42)
            energies2, configs2, accepted2 = nvt_monte_carlo(mc_routine, lattice, ham, 300.0, 50, 42)

            @test all(sum(c.components[1]) == sum(lattice.components[1]) for c in configs1)
            @test energies1 == energies2
            @test accepted1 == accepted2
        end
    end
end

@testset "nvt Monte Carlo atomistic version" begin
    @testset "nvt_monte_carlo_atomistic function tests" begin
        initial_config = FreeBirdIO.generate_random_starting_config(562.5, 2)
        at_init = AtomWalker(initial_config)

        lj = LJParameters(epsilon=0.1, sigma=2.5)

        # Metropolis Monte Carlo
        ats = LJAtomWalkers([at_init], lj)
        at = deepcopy(ats.walkers[1])

        cell_vec = at.configuration.cell.cell_vectors
        cell_volume = cell_vec[1][1] * cell_vec[2][2] * cell_vec[3][3]
        cell_size = cbrt(cell_volume).val

        # N = 4
        temperatures = collect(1000.0:-100:500.0)
        num_equilibration_steps = 1_000
        num_sampling_steps = 1_000
        step_size = cell_size / 50

        mc_params = MetropolisMCParameters(
            temperatures,
            equilibrium_steps=num_equilibration_steps,
            sampling_steps=num_sampling_steps,
            step_size=step_size,
            step_size_up=1.0,
            accept_range=(0.5,0.5)
        )

        mc_routine = MCRandomWalkMaxE()

        mc_energies, mc_ls, mc_cvs, acceptance_rates = monte_carlo_sampling(mc_routine, at, lj, mc_params)

        @test length(mc_energies) == length(temperatures)
        @test length(mc_ls.walkers) == length(temperatures)
        @test length(mc_cvs) == length(temperatures)
        @test length(acceptance_rates) == length(temperatures)

        @test all(isfinite, mc_energies)
        @test all(isfinite, mc_cvs)

        @test all(0 ≤ x ≤ 1 for x in acceptance_rates)

    end

    @testset "nvt Monte Carlo lattice version" begin
        lattice = SLattice{SquareLattice}(
            supercell_dimensions=(2, 2, 1),
            components=[[1,2]]
        )

        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

        mc_routine = MCNewSample()

        # Metropolis Monte Carlo
        energies, configs, accepted = nvt_monte_carlo(mc_routine, lattice, ham, 300.0, 100, 42)

        @test length(energies) == 100
        @test length(configs) == 100
        @test 0 ≤ accepted ≤ 100
        @test all(isfinite, energies)

        initial_lattice = SLattice{SquareLattice}(
            supercell_dimensions=(2, 2, 1),
            components=[[1,2]]
        )

        temperatures = collect(1000.0:-100:500.0)
        num_equilibration_steps = 1_000
        num_sampling_steps = 1_000

        mc_params = MetropolisMCParameters(
            temperatures,
            equilibrium_steps=num_equilibration_steps,
            sampling_steps=num_sampling_steps,
            step_size=0.1,
            step_size_up=0.2,
            accept_range=(0.3, 0.4)
        )

        mc_energies, mc_configs, mc_cvs, acceptance_rates = monte_carlo_sampling(mc_routine, initial_lattice, ham, mc_params)

        @test length(mc_energies) == length(temperatures)
        @test length(mc_configs) == length(temperatures)
        @test length(mc_cvs) == length(temperatures)
        @test length(acceptance_rates) == length(temperatures)

        @test all(isfinite, mc_energies)
        @test all(isfinite, mc_cvs)

        @test all(0 ≤ x ≤ 1 for x in acceptance_rates)

    end
end
@testset "NVT MC loop hardening (aliasing, null moves, phase seeding)" begin
    # Shared fixture: a 4-atom periodic LJ walker in a 10 A cubic box.
    # Calibration ledger: null-move surplus (accepted - changepoints)/n on
    #   MCMixedMoves(5, 1), n = 2000, T = 300 K, seeds 777/778/779: measured
    #   0.0435, 0.0380, 0.0435 (accepted identity swaps only; expectation
    #   (1/6)*(1/4) ~ 0.042, binomial sigma ~ 0.0045). The historical
    #   double-draw loop adds null moves at 5/36 ~ 0.139 on top, landing near
    #   0.146. Gate: surplus <= 0.09, over 10 sigma above the measured band's
    #   ceiling and over 12 sigma below the defect value.
    hardening_box = [[10.0u"Å", 0u"Å", 0u"Å"],
                     [0u"Å", 10.0u"Å", 0u"Å"],
                     [0u"Å", 0u"Å", 10.0u"Å"]]
    hardening_coords = [:H => [0.2, 0.2, 0.2], :H => [0.7, 0.3, 0.4],
                        :H => [0.3, 0.7, 0.6], :H => [0.6, 0.6, 0.3]]
    hardening_at = AtomWalker(FastSystem(periodic_system(hardening_coords, hardening_box, fractional=true)))
    hardening_lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.0, shift=true)

    # State-changing acceptances, counted from the stored snapshots
    function count_changepoints(configs, at_start)
        prev = at_start.configuration.position
        n = 0
        for c in configs
            if c.configuration.position != prev
                n += 1
            end
            prev = c.configuration.position
        end
        return n
    end

    @testset "Stored trajectory snapshots are independent (MCRandomWalkMaxE)" begin
        energies, configs, accepted = nvt_monte_carlo(
            MCRandomWalkMaxE(), deepcopy(hardening_at), hardening_lj, 300.0, 200, 0.4, 4242)
        @test 0 < accepted < 200
        # Distinct objects, not one aliased walker repeated
        @test all(configs[i] !== configs[j] for i in 1:20 for j in (i + 1):20)
        # The stored trajectory actually varies
        @test any(configs[i].configuration.position != configs[end].configuration.position for i in 1:199)
        # Each snapshot's energy field matches a from-scratch recompute of its
        # own configuration (incremental-accumulation drift class)
        for i in (1, 50, 100, 150, 200)
            E_re = interacting_energy(configs[i].configuration, hardening_lj,
                                      configs[i].list_num_par, configs[i].frozen) +
                   configs[i].energy_frozen_part
            @test isapprox(configs[i].energy, E_re; atol=1e-9u"eV")
            @test configs[i].energy == energies[i]
        end
        # Walk-only path: every acceptance changes the configuration, exactly
        @test accepted == count_changepoints(configs, hardening_at)
    end

    @testset "Null-move accounting (MCMixedMoves single channel draw)" begin
        n_steps = 2000
        energies, configs, accepted = nvt_monte_carlo(
            MCMixedMoves(5, 1), deepcopy(hardening_at), hardening_lj, 300.0, n_steps, 0.4, 777)
        # Accepted identity swaps (same atom drawn twice in the swap channel)
        # are the only legal accepted-without-change steps; the historical
        # double-draw loop added a 5/36 null-move channel on top. See the
        # calibration ledger above.
        surplus = (accepted - count_changepoints(configs, hardening_at)) / n_steps
        @test 0.0 <= surplus <= 0.09
        @test all(configs[i] !== configs[j] for i in 1:20 for j in (i + 1):20)
        for i in (1, 1000, 2000)
            @test configs[i].energy == energies[i]
        end
    end

    @testset "Equilibration and production streams are independent" begin
        # Lattice driver wiring: the production leg must reproduce an explicit
        # nvt_monte_carlo call seeded with random_seed + 1 from the
        # equilibration end configuration
        lattice = SLattice{SquareLattice}(supercell_dimensions=(2, 2, 1), components=[[1, 2]])
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        seed = 4242
        params = MetropolisMCParameters([300.0];
            equilibrium_steps=60, sampling_steps=80, random_seed=seed)
        d_energies, d_configs, d_cvs, d_rates = monte_carlo_sampling(MCNewSample(), lattice, ham, params)
        eq_e, eq_c, eq_acc = nvt_monte_carlo(MCNewSample(), lattice, ham, 300.0, 60, seed)
        pr_e, pr_c, pr_acc = nvt_monte_carlo(MCNewSample(), eq_c[end], ham, 300.0, 80, seed + 1)
        @test d_energies[1] == mean(pr_e)
        @test d_rates[1] == pr_acc / 80
        # Atomistic driver: same-seed reproducibility end to end
        p1 = MetropolisMCParameters([500.0];
            equilibrium_steps=100, sampling_steps=100, step_size=0.3, random_seed=99)
        p2 = MetropolisMCParameters([500.0];
            equilibrium_steps=100, sampling_steps=100, step_size=0.3, random_seed=99)
        e1, ls1, cv1, r1 = monte_carlo_sampling(MCRandomWalkMaxE(), deepcopy(hardening_at), hardening_lj, p1)
        e2, ls2, cv2, r2 = monte_carlo_sampling(MCRandomWalkMaxE(), deepcopy(hardening_at), hardening_lj, p2)
        @test e1 == e2
        @test cv1 == cv2
        @test r1 == r2
        @test ls1.walkers[1].configuration.position == ls2.walkers[1].configuration.position
    end
end
