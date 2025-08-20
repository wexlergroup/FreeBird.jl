@testset "AbstractLiveSets Tests" begin
    @testset "atomistic_livesets.jl tests" begin

        # Create periodic box and LJ parameters
        box = [[10.0u"Å", 0u"Å", 0u"Å"],
               [0u"Å", 10.0u"Å", 0u"Å"],
               [0u"Å", 0u"Å", 10.0u"Å"]]
        lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

        ljs = CompositeLJParameters(3, [lj for _ in 1:6]) # 3 components, 6 parameters

        # Test with 4 atoms: 2 H and 2 O atoms
        coor_list = [:H => [0.2, 0.5, 0.5],
                     :H => [0.4, 0.5, 0.5],
                     :O => [0.6, 0.5, 0.5],
                     :O => [0.8, 0.5, 0.5]]
        at = FastSystem(periodic_system(coor_list, box, fractional=true))

        surf_list = [:H => [0.0, 0.0, 0.0],
                     :H => [0.0, 0.5, 0.0],
                     :H => [0.5, 0.0, 0.0],
                     :H => [0.5, 0.5, 0.0]]

        surface = AtomWalker(FastSystem(periodic_system(surf_list, box, fractional=true)); freeze_species=[:H])
        surface.energy_frozen_part = interacting_energy(surface.configuration, lj)

        @testset "AtomWalker assign_energy funtion tests" begin
            
            # Test no frozen particles
            walker_free = AtomWalker(at)
            AbstractLiveSets.assign_energy!(walker_free, lj)
            @test walker_free.energy > 0.0u"eV"
            @test walker_free.energy_frozen_part == 0.0u"eV"
            
            # Test with O frozen
            walker_partial = AtomWalker(at, freeze_species=[:O])
            AbstractLiveSets.assign_frozen_energy!(walker_partial, lj)
            AbstractLiveSets.assign_energy!(walker_partial, lj)
            @test walker_partial.energy > walker_partial.energy_frozen_part > 0.0u"eV"
            
            # Test all frozen
            walker_frozen = AtomWalker(at, freeze_species=[:H, :O])
            AbstractLiveSets.assign_frozen_energy!(walker_frozen, lj)
            AbstractLiveSets.assign_energy!(walker_frozen, lj)
            @test walker_frozen.energy == walker_frozen.energy_frozen_part > 0.0u"eV"

            # Test with surface
            walker_surface = AtomWalker(at)
            AbstractLiveSets.assign_energy!(walker_surface, ljs, surface)
            @test walker_surface.energy ≠ 0.0u"eV"
            @test walker_surface.energy_frozen_part == 0.0u"eV"  # untouched
        end


        @testset "LJSurfaceWalkers struct and functions tests" begin
            
            # Create multiple walkers with different configurations
            walkers = [
                AtomWalker(at),  
                AtomWalker(at),  
                AtomWalker(at)
            ]

            @testset "Constructor with energy assignment" begin
                lj_walkers = LJSurfaceWalkers(walkers, ljs, surface; assign_energy=true)

                # Test type and structure
                @test lj_walkers isa LJSurfaceWalkers
                @test length(lj_walkers.walkers) == 3
                @test lj_walkers.potential === ljs
                @test lj_walkers.surface === surface
                
                # Test energy assignment
                @test all(w.energy > 0.0u"eV" for w in lj_walkers.walkers)
                @test lj_walkers.walkers[1].energy_frozen_part == interacting_energy(surface.configuration, lj)
                @test lj_walkers.walkers[1].energy_frozen_part == lj_walkers.walkers[2].energy_frozen_part == lj_walkers.walkers[3].energy_frozen_part  # All frozen particles
            
                # Test thread safety
                threaded_lj_walkers = LJSurfaceWalkers(walkers, ljs, surface, :threads)
                @test threaded_lj_walkers isa LJSurfaceWalkers
                @test length(threaded_lj_walkers.walkers) == 3
                @test threaded_lj_walkers.potential === ljs
                @test threaded_lj_walkers.surface === surface
                @test all(threaded_lj_walkers.walkers[i].energy == lj_walkers.walkers[i].energy for i in 1:3)
                @test all(threaded_lj_walkers.walkers[i].energy_frozen_part == lj_walkers.walkers[i].energy_frozen_part for i in 1:3)

                # Test distributed safety
                distributed_lj_walkers = LJSurfaceWalkers(walkers, ljs, surface, :distributed)
                @test distributed_lj_walkers isa LJSurfaceWalkers
                @test length(distributed_lj_walkers.walkers) == 3
                @test distributed_lj_walkers.potential === ljs
                @test distributed_lj_walkers.surface === surface
                @test all(distributed_lj_walkers.walkers[i].energy == lj_walkers.walkers[i].energy for i in 1:3)
                @test all(distributed_lj_walkers.walkers[i].energy_frozen_part == lj_walkers.walkers[i].energy_frozen_part for i in 1:3)
            end

            @testset "Constructor without energy assignment" begin
                # Reset walker energies
                for w in walkers
                    w.energy = 0.0u"eV"
                    w.energy_frozen_part = 0.0u"eV"
                end

                lj_walkers = LJSurfaceWalkers(walkers, ljs, surface, assign_energy=false)

                # Verify energies weren't assigned
                @test all(w.energy == 0.0u"eV" for w in lj_walkers.walkers)
                @test all(w.energy_frozen_part == 0.0u"eV" for w in lj_walkers.walkers)
            end
        end

        @testset "LJAtomWalkers struct and functions tests" begin
            
            # Create multiple walkers with different configurations
            walkers = [
                AtomWalker(at),  # No frozen atoms
                AtomWalker(at, freeze_species=[:O]),  # O atoms frozen
                AtomWalker(at, freeze_species=[:H, :O])  # All frozen
            ]

            @testset "Constructor with energy assignment" begin
                lj_walkers = LJAtomWalkers(walkers, lj; const_frozen_part=false)
                
                # Test type and structure
                @test lj_walkers isa LJAtomWalkers
                @test length(lj_walkers.walkers) == 3
                @test lj_walkers.potential === lj
                
                # Test energy assignment
                @test all(w.energy > 0.0u"eV" for w in lj_walkers.walkers)
                @test lj_walkers.walkers[1].energy_frozen_part == 0.0u"eV"  # No frozen particles
                @test lj_walkers.walkers[2].energy_frozen_part > 0.0u"eV"   # Some frozen
                @test lj_walkers.walkers[3].energy == lj_walkers.walkers[3].energy_frozen_part  # All frozen
            end

            @testset "Constructor without energy assignment" begin
                # Reset walker energies
                for w in walkers
                    w.energy = 0.0u"eV"
                    w.energy_frozen_part = 0.0u"eV"
                end
                
                lj_walkers = LJAtomWalkers(walkers, lj, assign_energy=false)
                
                # Verify energies weren't assigned
                @test all(w.energy == 0.0u"eV" for w in lj_walkers.walkers)
                @test all(w.energy_frozen_part == 0.0u"eV" for w in lj_walkers.walkers)
            end

            @testset "LJAtomWalkers show method tests" begin
                # Use existing setup from parent testset
                walkers = [
                    AtomWalker(at),
                    AtomWalker(at, freeze_species=[:O]),
                    AtomWalker(at, freeze_species=[:H, :O])
                ]
                
                # Test display with small number of walkers
                lj_walkers_small = LJAtomWalkers(walkers, lj)
                output_small = sprint(show, lj_walkers_small)
                
                # Verify key components are shown
                @test occursin("LJAtomWalkers", output_small)
                @test all(occursin("[$i]", output_small) for i in 1:3)
                @test occursin("LJParameters", output_small)
                
                # Test display with truncation (>10 walkers)
                walkers_large = [AtomWalker(at) for _ in 1:15]
                lj_walkers_large = LJAtomWalkers(walkers_large, lj)
                output_large = sprint(show, lj_walkers_large)
                
                @test occursin("Omitted 5 walkers", output_large)
                @test occursin("[1]", output_large)
                @test occursin("[15]", output_large)
            end

            @testset "LJSurfaceWalkers show method tests" begin
                # Use existing setup from parent testset
                walkers = [
                    AtomWalker(at),
                    AtomWalker(at, freeze_species=[:O]),
                    AtomWalker(at, freeze_species=[:H, :O])
                ]
                
                # Test display with small number of walkers
                lj_walkers_small = LJSurfaceWalkers(walkers, ljs, surface; assign_energy=true)
                output_small = sprint(show, lj_walkers_small)
                
                # Verify key components are shown
                @test occursin("LJSurfaceWalkers", output_small)
                @test all(occursin("[$i]", output_small) for i in 1:3)
                @test occursin("LJParameters", output_small)
                @test occursin("Surface: ", output_small)
                
                # Test display with truncation (>10 walkers)
                walkers_large = [AtomWalker(at) for _ in 1:15]
                lj_walkers_large = LJSurfaceWalkers(walkers_large, ljs, surface; assign_energy=true)
                output_large = sprint(show, lj_walkers_large)
                
                @test occursin("Omitted 5 walkers", output_large)
                @test occursin("[1]", output_large)
                @test occursin("[15]", output_large)
            end
        end
    end


    @testset "lattice_livesets.jl tests" begin
        @testset "LatticeWalker assign_energy function tests" begin
            # Setup a simple lattice configuration for testing
            square_lattice = MLattice{1,SquareLattice}(
                        lattice_constant=1.0,                   # Unit cell size
                        basis=[(0.0, 0.0, 0.0)],                # Single atom per unit cell
                        supercell_dimensions=(4, 4, 1),         # 4x4x1 supercell
                        periodicity=(true, true, false),        # Periodic in x,y but not z
                        cutoff_radii=[1.1, 1.5],                # Neighbor cutoff distances
                        components=:equal,                      # Split into equal components
                        adsorptions=:full                       # All sites are adsorption sites
                    )
            
            # Create a simple classical hamiltonian for testing
            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
            
            @testset "Basic energy assignment" begin
                # Test without energy perturbation
                s_walker = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )
                AbstractLiveSets.assign_energy!(s_walker, ham)
                
                # Test that energy has been assigned with correct units
                @test !iszero(s_walker.energy)
                @test unit(s_walker.energy) == u"eV"
                
                # Store initial energy for comparison
                initial_energy = s_walker.energy
                
                # Test that repeated assignments with same configuration give same energy
                AbstractLiveSets.assign_energy!(s_walker, ham)
                @test s_walker.energy == initial_energy
            end
        
            @testset "Energy perturbation" begin
                # Test with energy perturbation
                s_walker1 = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )
                s_walker2 = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )
                
                # Assign energy with perturbation
                perturb_amount = 1.0
                AbstractLiveSets.assign_energy!(s_walker1, ham, perturb_energy=perturb_amount)
                AbstractLiveSets.assign_energy!(s_walker2, ham, perturb_energy=perturb_amount)
                
                # Test that perturbed energies are different due to random perturbation
                @test s_walker1.energy ≠ s_walker2.energy
                
                # Test that perturbation is within expected bounds
                base_energy = interacting_energy(square_lattice, ham)
                max_perturbation = perturb_amount * 0.5u"eV"  # Based on implementation using (rand() - 0.5)
                
                @test abs(s_walker1.energy - base_energy) ≤ max_perturbation
                @test abs(s_walker2.energy - base_energy) ≤ max_perturbation
            end
        
            @testset "Zero perturbation" begin
                # Test that zero perturbation gives consistent results
                s_walker = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )
                AbstractLiveSets.assign_energy!(s_walker, ham, perturb_energy=0.0)
                
                expected_energy = interacting_energy(square_lattice, ham)
                @test s_walker.energy == expected_energy
            end
            
            @testset "Type stability" begin
                s_walker = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )
                # Test that energy always maintains eV units
                AbstractLiveSets.assign_energy!(s_walker, ham, perturb_energy=1.0)
                @test typeof(s_walker.energy) == typeof(1.0u"eV")
                
                AbstractLiveSets.assign_energy!(s_walker, ham, perturb_energy=0.0)
                @test typeof(s_walker.energy) == typeof(1.0u"eV")
            end
        end


        @testset "LatticeGasWalkers struct and functions tests" begin
            
            # Setup test components
            lattice = MLattice{2,SquareLattice}(
                lattice_constant=1.0,
                supercell_dimensions=(2, 2, 1),
                components=:equal,
                cutoff_radii=[1.1, 1.5]
            )

            hams_2x2 = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4]
            mlham = MLatticeHamiltonian(2, hams_2x2)
            walkers = [LatticeWalker(lattice) for _ in 1:3]
        
            @testset "Constructor and Energy Assignment" begin
                # Test without energy assignment
                lgw_no_energy = LatticeGasWalkers(walkers, mlham, assign_energy=false)
                @test all(w.energy == 0.0u"eV" for w in lgw_no_energy.walkers)
                
                # Test with energy assignment (no perturbation)
                lgw = LatticeGasWalkers(walkers, mlham, assign_energy=true)
                @test all(w.energy != 0.0u"eV" for w in lgw.walkers)
                @test all(typeof(w.energy) == typeof(1.0u"eV") for w in lgw.walkers)
                
                # Test energy assignment with perturbation
                perturb_amount = 1.0
                lgw_perturbed = LatticeGasWalkers(walkers, mlham, assign_energy=true, perturb_energy=perturb_amount)
                perturbed_energies = [w.energy for w in lgw_perturbed.walkers]
                
                # Check energies are different when perturbed
                @test length(unique(perturbed_energies)) == length(perturbed_energies)
                
                # Check perturbation bounds
                base_energy = interacting_energy(lattice, mlham)
                @test all(abs(e - base_energy) ≤ perturb_amount * 0.5u"eV" for e in perturbed_energies)
            end

            @testset "Show method tests" begin
                # Test both small and large sets of walkers
                small_walkers = [LatticeWalker(lattice) for _ in 1:3]
                large_walkers = [LatticeWalker(lattice) for _ in 1:15]

                @test occursin("Omitted", sprint(show, large_walkers))
                
                lgw_small = LatticeGasWalkers(small_walkers, mlham)
                lgw_large = LatticeGasWalkers(large_walkers, mlham)
                
                small_output = sprint(show, lgw_small)
                large_output = sprint(show, lgw_large)

                # Check key components
                @test occursin("LatticeGasWalkers", small_output)
                @test occursin("lattice_vectors", small_output)
                @test occursin("supercell_dimensions", small_output)
                @test occursin("periodicity", small_output)
                @test occursin("basis", small_output)
                
                @test all(occursin("[$i]", small_output) for i in 1:3)

                # Check truncation for large set
                @test occursin("Omitted 5 walkers", large_output)
                @test occursin("energy =", large_output)
                @test occursin("iter =", large_output)
                @test occursin("[15]", large_output)
            end
        end


        @testset "print_lattice_walker_in_walkers function tests" begin
            # Setup a simple lattice configuration for testing
            square_lattice = MLattice{2,SquareLattice}(
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

            # Test printing output
            output = sprint(io -> AbstractLiveSets.print_lattice_walker_in_walkers(io, s_walker))

            # Check output format for MLattice
            @test occursin("occupations:", output)
            @test occursin("components:", output)
            @test occursin("component 1:", output)
            @test occursin("component 2:", output)
            
            # Test component values are printed
            for (i, comp) in enumerate(s_walker.configuration.components)
                @test occursin("$comp", output)
            end
        end
    end

end