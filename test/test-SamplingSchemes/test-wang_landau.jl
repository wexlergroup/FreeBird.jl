@testset "wang_landau.jl tests" begin
    @testset "wang_landau function tests" begin

        # Setup
        lattice = SLattice{SquareLattice}(
            supercell_dimensions=(4, 4, 1),
            components=[[1]]
        )

        
        adsorption_energy = -0.04
        nn_energy = -0.01
        nnn_energy = -0.0025

        energy_min = adsorption_energy - 0.01
        energy_max = adsorption_energy + 0.01
        num_energy_bins = 10

        ham = GenericLatticeHamiltonian(adsorption_energy, [nn_energy, nnn_energy], u"eV")

        wl_params = WangLandauParameters(energy_min=energy_min, energy_max=energy_max, num_energy_bins=num_energy_bins)

        @testset "Basic functionality" begin
            energies_wl, configs, wl_params, S, H = wang_landau(
                lattice, ham, wl_params
            )
            
            @test length(S) == num_energy_bins
            @test length(H) == num_energy_bins
            @test length(energies_wl) == 28
            @test length(configs) == 28
            @test all(x -> x ≥ 0, H)
        end

        @testset "Energy conservation" begin
            energy_bins = range(-0.1, 0.1, length=20)
            wl_params = WangLandauParameters(50, 0.8, 2.7, 1.001, energy_bins, 1000, 0.01, 42)
            energies_wl, configs, wl_params, S, H = wang_landau(
                lattice, ham, wl_params
            )
            
            @test all(e -> minimum(energy_bins) ≤ e ≤ maximum(energy_bins), energies_wl)
            @test all(c -> sum(c.components[1]) == sum(lattice.components[1]), configs)
        end

        @testset "Modification factor" begin
            energy_bins = range(-0.1, 0.1, length=20)
            wl_params = WangLandauParameters(50, 0.8, 2.7, 1.1, energy_bins, 1000, 0.01, 42)
            _, _, _, S1, _ = wang_landau(
                lattice, ham, wl_params
            )
            wl_params = WangLandauParameters(50, 0.8, 2.7, 1.001, energy_bins, 1000, 0.01, 42)
            _, _, _, S2, _ = wang_landau(
                lattice, ham, wl_params
            )
            
            @test !isapprox(S1, S2)
        end
    end

    @testset "Wang Landau atomistic sampling" begin
        # Setup test system
        box = [[10.0u"Å", 0u"Å", 0u"Å"],
                    [0u"Å", 10.0u"Å", 0u"Å"],
                    [0u"Å", 0u"Å", 10.0u"Å"]]
        lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
        
        coor_list = [:H => [0.2, 0.5, 0.5],
                        :H => [0.4, 0.5, 0.5]]
        at = FastSystem(periodic_system(coor_list, box, fractional=true))
        
        walkers = [AtomWalker(at) for _ in 1:3]
        liveset = LJAtomWalkers(walkers, lj)

        wk = liveset.walkers[1]

        wl_params = WangLandauParameters(f_min=1.01, energy_min=-1.27, energy_max=10.0, num_energy_bins=1000, max_iter=1)

        wl_energies, wl_configs, wl_params, S, H = wang_landau(wk, lj, wl_params)

        @test wl_energies isa Vector{Float64}
        @test wl_configs isa Vector{AtomWalker}
        @test length(S) == 1000
        @test length(H) == 1000

        wl_params = WangLandauParameters(f_min=1.01, energy_min=-1.27, energy_max=15.0, num_energy_bins=2)
        wl_energies, wl_configs, wl_params, S, H = wang_landau(wk, lj, wl_params)

        @test wl_energies isa Vector{Float64}
        @test wl_configs isa Vector{AtomWalker}
        @test length(S) == 2
        @test length(H) == 2

        ls = LJAtomWalkers([AtomWalker{1}(config.configuration) for config in wl_configs], lj)


    end


    @testset "get_bin_index function tests" begin
        bins = [-1.0, -0.5, 0.0, 0.5, 1.0]

        @test SamplingSchemes.get_bin_index(-0.7, bins) == 1
        @test SamplingSchemes.get_bin_index(-0.2, bins) == 2
        @test SamplingSchemes.get_bin_index(0.3, bins) == 3
        @test SamplingSchemes.get_bin_index(0.8, bins) == 4
        @test SamplingSchemes.get_bin_index(1.5, bins) == 5  # Edge case: above max, insert at end
        @test SamplingSchemes.get_bin_index(0.0, bins) == 2  # Edge case: exact boundary, put in bin 2
    end

    @testset "AtomicLattice lower energy boundary" begin
        lattice = AtomicLattice{1,SquareLattice}(
            lattice_atom="Pd", surface=:fcc100,
            supercell_dimensions=(2, 2, 1), lattice_constant=3.947,
            periodicity=(true, true, false), adsorbate_atoms=["O"],
            components=[1], num_nearest_neighbors=0,
            type_of_sites=["hollow"])
        h = GenericLatticeHamiltonian(-0.04, Float64[], u"eV")
        @test SamplingSchemes._lattice_bin_index(lattice, -0.04, [-0.04, 0.0]) == 1
        @test SamplingSchemes.get_bin_index(-0.04, [-0.04, 0.0]) == 0
        parameters = WangLandauParameters(
            num_steps=2, f_min=2.0, energy_min=-0.04,
            energy_max=-0.03, num_energy_bins=3, max_iter=0,
            random_seed=42)
        energies, configurations, _, entropy, histogram =
            wang_landau(lattice, h, parameters)
        @test first(energies) == -0.04
        @test configurations[1] isa AtomicLattice
        @test entropy[1] > 0
        @test all(>=(0), histogram)

        for bins in (Float64[], [-0.04], [-0.04, -0.04],
                     [-0.03, -0.04], [NaN, 0.0], [-0.05, Inf])
            invalid = WangLandauParameters(1, 0.8, 2.0, 1.5, bins, 0, 0.01, 42)
            @test_throws ArgumentError wang_landau(lattice, h, invalid)
        end
        for (field, values) in ((:num_steps, (0, -1)), (:max_iter, (-1,)),
                               (:flatness_criterion, (-0.1, 1.0, NaN, Inf)),
                               (:f_initial, (0.0, 1.0, NaN, Inf)),
                               (:f_min, (0.0, 1.0, NaN, Inf)))
            for value in values
                invalid = deepcopy(parameters)
                setproperty!(invalid, field, value)
                @test_throws ArgumentError wang_landau(lattice, h, invalid)
            end
        end
        outside_window = WangLandauParameters(
            num_steps=1, f_min=2.0, energy_min=0.0,
            energy_max=1.0, num_energy_bins=2, max_iter=0)
        @test_throws ArgumentError wang_landau(lattice, h, outside_window)
        namespace = FreeBird.EnergyEval.pydict()
        FreeBird.EnergyEval.pyimport("builtins").exec("""
class FreeBirdNonfiniteWLCalculator:
    def get_potential_energy(self, atoms): return float('nan')
""", namespace)
        nonfinite = PyMLPotential(FreeBird.AbstractPotentials.ASEcalculator(
            namespace["FreeBirdNonfiniteWLCalculator"]()))
        @test_throws DomainError wang_landau(lattice, nonfinite, parameters)
    end

    if isfile("entropy.csv")
        rm("entropy.csv", force=true)
    end
end
