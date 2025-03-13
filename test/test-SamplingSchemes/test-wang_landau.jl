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
        at = generate_initial_configs(1, 562.5, 2; particle_type=:H)[1]
        wk = AtomWalker(at)
        lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

        wl_params = WangLandauParameters(f_min=1.01, energy_min=-1.27, energy_max=10.0, num_energy_bins=1000)

        wl_energies, wl_configs, wl_params, S, H = wang_landau(wk, lj, wl_params)

        ls = LJAtomWalkers([AtomWalker{1}(config.configuration) for config in wl_configs], lj)

        @test wl_energies isa Vector{Float64}
        @test wl_configs isa Vector{AtomWalker}
        @test length(S) == 1000
        @test length(H) == 1000
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
end