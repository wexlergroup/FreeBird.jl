using Test
using StaticArrays
using Random
using Unitful

@testset "Monte Carlo Moves tests" begin
    
    @testset "random_walks.jl tests" begin

        @testset "single_atom_random_walk function" begin

            # Test basic functionality
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 2.0, 3.0] .* u"Å")
            step_size = 0.5
            new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size)
            
            @test new_pos ≠ pos
            
            # Test step size constraints
            for _ in 1:10
                pos = SVector{3, typeof(0.0u"Å")}([0.0, 0.0, 0.0] .* u"Å")
                new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size)
                
                # Check that no component moved more than step_size
                @test all(abs.(new_pos - pos) .<= step_size * u"Å")
            end
            
            # Test type stability
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 1.0, 1.0] .* u"Å")
            @test typeof(MonteCarloMoves.single_atom_random_walk!(pos, step_size)) == typeof(pos)

        end

        @testset "MC_random_walk functions" begin

            @testset "AtomWalker MC_random_walk" begin

                # Setup test parameters
                n_steps = 100
                step_size = 0.5
                emax = 1.0u"eV"
                
                # Create a simple test system
                at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[2,1,3,4,5,6];particle_types=[:H,:O,:H,:Fe,:Au,:Cl])    # Example values
                atw = AtomWalker(at;freeze_species=[:H],merge_same_species=true)
                lj = LJParameters()  # Create test LJ parameters
                
                # Run MC simulation
                accept_this_walker, accept_rate, updated_atw = MC_random_walk!(
                    n_steps, atw, lj, step_size, emax
                )
                
                # Tests
                @test typeof(accept_this_walker) == Bool
                @test 0.0 <= accept_rate <= 1.0
                @test typeof(updated_atw) == typeof(atw)
                @test updated_at.energy <= emax
            end

            @testset "LatticeWalker MC_random_walk" begin

                # Setup test parameters
                n_steps = 100
                step_size = 0.5
                emax = 1.0
                
                # Create a simple square lattice system for test
                square_lattice = MLattice{1,SquareLattice}(
                    lattice_constant=1.0,                   # Unit cell size
                    basis=[(0.0, 0.0, 0.0)],                # Single atom per unit cell
                    supercell_dimensions=(4, 4, 1),         # 4x4x1 supercell
                    periodicity=(true, true, false),        # Periodic in x,y but not z
                    cutoff_radii=[1.1, 1.5],                # Neighbor cutoff distances
                    components=:equal,                      # Split into equal components
                    adsorptions=:full                       # All sites are adsorption sites
                )

                s_walker = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )

                # Create a simple triangular lattice system for test
                triangular_lattice = MLattice{1,TriangularLattice}(
                    lattice_constant=1.0,
                    basis=[(0.0, 0.0, 0.0), (1/2, sqrt(3)/2, 0.0)],  # Two atoms per unit cell
                    supercell_dimensions=(4, 2, 1),
                    periodicity=(true, true, false),
                    cutoff_radii=[1.1, 1.5],
                    components=:equal,
                    adsorptions=:full
                )

                t_walker = LatticeWalker(
                    triangular_lattice,
                    energy=5.0u"eV",
                    iter=0
                )

                # Create a simple Hamiltonian for testing
                h = LatticeGasHamiltonian(1u"eV", 2.2u"eV",3.33u"eV")

                # Test basic functionality
                @testset "Basic MC walk" begin

                    # Store initial state
                    initial_config = deepcopy(s_walker.configuration)
                    initial_energy = s_walker.energy
                    
                    # Perform MC walk
                    accepted, rate, updated_walker = MC_random_walk!(n_steps, s_walker, h, emax)            # !!!Error: type MLattice has no field occupations
                    
                    # Basic checks
                    @test typeof(accepted) == Bool
                    @test 0.0 <= rate <= 1.0
                    @test updated_walker isa LatticeWalker
                end

            end

        end

        @testset "MC_new_sample function" begin

            # Create a simple square lattice system for test
            square_lattice = MLattice{1,SquareLattice}(
                lattice_constant=1.0,                   # Unit cell size
                basis=[(0.0, 0.0, 0.0)],                # Single atom per unit cell
                supercell_dimensions=(4, 4, 1),         # 4x4x1 supercell
                periodicity=(true, true, false),        # Periodic in x,y but not z
                cutoff_radii=[1.1, 1.5],                # Neighbor cutoff distances
                components=:equal,                      # Split into equal components
                adsorptions=:full                       # All sites are adsorption sites
            )

            s_walker = LatticeWalker(
                square_lattice,
                energy=5.0u"eV",
                iter=0
            )

            

        end


    end

    @testset "nve_walks.jl tests" begin

    end

    @testset "helpers.jl tests" begin

    end

end

