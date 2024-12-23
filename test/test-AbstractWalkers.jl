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



end

