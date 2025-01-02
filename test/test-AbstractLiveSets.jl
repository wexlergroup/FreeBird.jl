@testset "AbstractLiveSets Tests" begin
    @testset "atomistic_livesets.jl tests" begin

        # Create periodic box and LJ parameters
        box = [[10.0u"Å", 0u"Å", 0u"Å"],
               [0u"Å", 10.0u"Å", 0u"Å"],
               [0u"Å", 0u"Å", 10.0u"Å"]]
        lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

        # Test with 4 atoms: 2 H and 2 O atoms
        coor_list = [:H => [0.2, 0.5, 0.5],
                     :H => [0.4, 0.5, 0.5],
                     :O => [0.6, 0.5, 0.5],
                     :O => [0.8, 0.5, 0.5]]
        at = FastSystem(periodic_system(coor_list, box, fractional=true))


        @testset "assign_energy funtion tests" begin
            
            # Test no frozen particles
            walker_free = AtomWalker(at)
            AbstractLiveSets.assign_energy!(walker_free, lj)
            @test walker_free.energy > 0.0u"eV"
            @test walker_free.energy_frozen_part == 0.0u"eV"
            
            # Test with O frozen
            walker_partial = AtomWalker(at, freeze_species=[:O])
            AbstractLiveSets.assign_energy!(walker_partial, lj)
            @test walker_partial.energy > walker_partial.energy_frozen_part > 0.0u"eV"
            
            # Test all frozen
            walker_frozen = AtomWalker(at, freeze_species=[:H, :O])
            AbstractLiveSets.assign_energy!(walker_frozen, lj)
            @test walker_frozen.energy == walker_frozen.energy_frozen_part > 0.0u"eV"
        end


        @testset "LJAtomWalkers struct and functions tests" begin
            
            # Create multiple walkers with different configurations
            walkers = [
                AtomWalker(at),  # No frozen atoms
                AtomWalker(at, freeze_species=[:O]),  # O atoms frozen
                AtomWalker(at, freeze_species=[:H, :O])  # All frozen
            ]

            @testset "Constructor with energy assignment" begin
                lj_walkers = LJAtomWalkers(walkers, lj)
                
                # Test type and structure
                @test lj_walkers isa LJAtomWalkers
                @test length(lj_walkers.walkers) == 3
                @test lj_walkers.lj_potential === lj
                
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
        end
    end

end