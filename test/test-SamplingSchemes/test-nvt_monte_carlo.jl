@testset "nvt_monte_carlo.jl tests" begin
    @testset "nvt_monte_carlo function tests" begin
        
        # Setup
        lattice = SLattice{SquareLattice}(
                supercell_dimensions=(2, 1, 1),
                components=[[true, false]]
            )
        
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

        @testset "Basic functionality" begin
            
            energies, configs, accepted = nvt_monte_carlo(lattice, ham, 300.0, 100, 42)
            
            @test length(energies) == 100
            @test length(configs) == 100
            @test 0 ≤ accepted ≤ 100
            @test all(isfinite, energies)
        end
    
        @testset "Temperature effects" begin
    
            energies_cold, _, accepted_cold = nvt_monte_carlo(lattice, ham, 1.0, 50, 42)
            energies_hot, _, accepted_hot = nvt_monte_carlo(lattice, ham, 1000.0, 50, 42)
    
            @test std(energies_cold) ≤ std(energies_hot)
            @test accepted_cold ≤ accepted_hot
        end
    
        @testset "Conservation & Reproducibility" begin
    
            energies1, configs1, accepted1 = nvt_monte_carlo(lattice, ham, 300.0, 50, 42)
            energies2, configs2, accepted2 = nvt_monte_carlo(lattice, ham, 300.0, 50, 42)
    
            @test all(sum(c.components[1]) == sum(lattice.components[1]) for c in configs1)
            @test energies1 == energies2
            @test accepted1 == accepted2
        end
    end
end