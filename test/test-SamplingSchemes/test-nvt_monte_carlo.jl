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