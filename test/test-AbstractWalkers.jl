@testset "AbstractWalkers tests" begin
    
    @testset "AtomWalkers.jl test" begin

        @testset "AtomWalker structure" begin
            at = FreeBirdIO.generate_multi_type_random_starting_config(10.0, [2,1], particle_types=[:H,:O])

            # Valid construction
            walker = AtomWalker{2}(at, list_num_par=[2,1])
            @test walker isa AtomWalker{2}
            @test walker.configuration == at
            @test walker.energy == 0.0u"eV"
            @test walker.iter == 0
            @test walker.frozen == [false, false]
            @test walker.energy_frozen_part == 0.0u"eV"

            # Custom parameters
            walker = AtomWalker{2}(at, 
                energy=1.0u"eV", 
                iter=5, 
                list_num_par=[2,1], 
                frozen=[true,false]
            )
            @test walker.energy == 1.0u"eV"
            @test walker.iter == 5
            @test walker.list_num_par == [2,1]
            @test walker.frozen == [true,false]        

            # Error outputs
            @test_throws ArgumentError AtomWalker{3}(at)                        # Wrong number of components
            @test_throws ArgumentError AtomWalker{2}(at, list_num_par=[1])      # Incorrect list_num_par length
            @test_throws ArgumentError AtomWalker{2}(at, frozen=[true])         # Incorrect frozen length
            @test_throws ArgumentError AtomWalker{2}(at, list_num_par=[1,1])    # Incorrect sum of list_num_par

            # Vector{AtomWalker}
            walkers = [AtomWalker(at), AtomWalker(at)]
            output = sprint(show, walkers)
            @test occursin("Vector{AtomWalker{2}}(2):", output)
            @test occursin("[1]", output)
            @test occursin("[2]", output)

        end
        
        
        @testset "AtomWalker function" begin
            # Create an example system for test
            at = FreeBirdIO.generate_multi_type_random_starting_config(
                10.0, [2,1,3,4,5,6], 
                particle_types=[:H,:O,:H,:Fe,:Au,:Cl]
            )
            
            # Use default parameters
            walker = AtomWalker(at)
            @test walker isa AtomWalker{5}  # Merged same species
            @test all(.!walker.frozen)      # No frozen components

            # Check components sorting
            components = split_components(walker.configuration, walker.list_num_par)
            atomic_numbers = [atomic_number(c)[1] for c in components]
            @test issorted(atomic_numbers)

            # Without merging same species
            walker_unmerged = AtomWalker(at, merge_same_species=false)
            @test walker_unmerged isa AtomWalker{6}  # Distinct components 

            # Freeze single species
            walker_frozen1 = AtomWalker(at, freeze_species=[:H])
            @test any(walker_frozen1.frozen)
            @test count(walker_frozen1.frozen) == 1

            # Freeze multiple species
            walker_frozen2 = AtomWalker(at, freeze_species=[:H, :O])
            @test any(walker_frozen2.frozen)
            @test count(walker_frozen2.frozen) == 2

            # Invalid freeze species
            @test_throws ArgumentError AtomWalker(at, freeze_species=[:X])

        end

        @testset "update_walker function" begin
            at = FreeBirdIO.generate_multi_type_random_starting_config(
                10.0, [2,1,4,5,6], 
                particle_types=[:H,:O,:Fe,:Au,:Cl]
            )

            walker = AtomWalker(at)
            count_frozen = count(walker.frozen)
            
            # Test updating energy
            new_energy = 1.5u"eV"
            updated = update_walker!(walker, :energy, new_energy)
            @test updated.energy == new_energy
            @test walker.energy == new_energy  # Test original walker modified
            
            # Test updating iteration count
            update_walker!(walker, :iter, 10)
            @test walker.iter == 10
            
            # Test updating frozen status
            new_frozen = [true, false, false, false, false,]
            update_walker!(walker, :frozen, new_frozen)
            @test walker.frozen == new_frozen
            @test count(walker.frozen) != count_frozen

            # Test invalid property
            @test_throws ErrorException update_walker!(walker, :invalid_prop, 0)
        
        end
    end

    @testset "LatticeWalkers.jl test" begin

        @testset "compute_neighbors function" begin
            # Simple cubic lattice
            lattice = [
                2.0 0.0 0.0;
                0.0 2.0 0.0;
                0.0 0.0 2.0
            ]

            positions = [
                0.0 0.0 0.0;
                1.0 1.0 0.0;
                0.0 1.0 1.0;
                1.0 0.0 1.0
                1.0 1.0 1.0
            ]

            periodicity = (false, false, false)

            cutoff_radii = [1.0, 1.5]

            neighbors = compute_neighbors(lattice, positions, periodicity, cutoff_radii)

            @test sort(neighbors[1][1]) == []
            @test sort(neighbors[1][2]) == [2, 3, 4]
            @test sort(neighbors[2][1]) == [5]
            @test sort(neighbors[2][2]) == [1, 3, 4]            
            @test sort(neighbors[5][1]) == [2, 3, 4]
            @test sort(neighbors[5][2]) == []

            # Periodic boundary conditions
            lattice = [
                2.0 0.0 0.0;
                0.0 2.0 0.0;
                0.0 0.0 2.0
            ]

            positions = [
                0.0 0.0 0.0;
                1.9 0.0 0.0   
            ]

            periodicity = (true, false, false)
            cutoff_radii = [1.5, 2.0]

            neighbors = compute_neighbors(lattice, positions, periodicity, cutoff_radii)

            @test neighbors[1][1] == [2]
            @test length(neighbors[1][2]) == 0

        end

        @testset "lattice_positions function" begin
            # Simple cubic lattice
            lattice_vectors = [
                1.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0
            ]

            basis = [(0.0, 0.0, 0.0)]
            dims = (2, 2, 1)

            positions = lattice_positions(lattice_vectors, basis, dims)

            @test size(positions) == (4, 3)
            @test positions[1,:] ≈ [0.0, 0.0, 0.0]
            @test positions[2,:] ≈ [1.0, 0.0, 0.0]
            @test positions[3,:] ≈ [0.0, 1.0, 0.0]
            @test positions[4,:] ≈ [1.0, 1.0, 0.0]

            # FCC lattice
            lattice_vectors = [
                1.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0
            ]
            basis = [
                (0.0, 0.0, 0.0),
                (0.5, 0.5, 0.0),
                (0.5, 0.0, 0.5),
                (0.0, 0.5, 0.5)
            ]
            dims = (1, 1, 1)
            
            positions = lattice_positions(lattice_vectors, basis, dims)
            
            @test size(positions) == (4, 3)
            for pos in eachrow(positions)
                @test 0.0 <= pos[1] <= 1.0
                @test 0.0 <= pos[2] <= 1.0
                @test 0.0 <= pos[3] <= 1.0
            end
        
            # Non-orthogonal lattice
                lattice_vectors = [
                    1.0 0.5 0.0
                    0.0 1.0 0.0
                    0.0 0.0 1.0
                ]

                basis = [(0.0, 0.0, 0.0)]
                dims = (2, 2, 1)
                
                positions = lattice_positions(lattice_vectors, basis, dims)
                
                @test size(positions) == (4, 3)
                @test positions[2,1] > positions[1,1]
                @test positions[3,2] > positions[1,2]
            

        end
        
        @testset "split_into_subarrays function" begin
            # 
        end

        @testset "mlattice_setup function" begin
            #
        end

        @testset "MLattice function" begin
            # Test default construction
            lattice = MLattice{2,SquareLattice}()
            @test size(lattice.lattice_vectors) == (3,3)
            @test lattice.supercell_dimensions == (4,4,1)
            @test length(lattice.components) == 2

            # Test custom parameters
            lattice = MLattice{3,SquareLattice}(
                lattice_constant=2.0,
                supercell_dimensions=(6,6,1),
                components=[[1],[2],[3,4]],
                adsorptions=[1,2]
            )

            @test lattice.lattice_vectors[1,1] == 2.0
            @test length(lattice.components) == 3
            @test count(lattice.adsorptions) == 2

            # Test triangular lattice
            lattice = MLattice{2,TriangularLattice}()
            @test size(lattice.lattice_vectors) == (3,3)
            @test length(lattice.basis) == 2
            @test lattice.supercell_dimensions == (4,2,1)
            
            @test_throws BoundsError MLattice{2,TriangularLattice}(
                components=[[1,2]], # Wrong number of components
                adsorptions=:full
            )
        end

        @testset "LatticeWalker structure" begin
            #
        end

    end

end
