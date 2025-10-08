@testset "Lattice live sets tests" begin
    @testset "LatticeGasWalkers with SLattice tests" begin
        # Setup a simple lattice configuration for testing
        square_lattice = SLattice{SquareLattice}(supercell_dimensions=(4, 4, 1), components=[[1,2,3,4]]) # 4x4 square lattice with one component

        lattices = [generate_random_new_lattice_sample!(deepcopy(square_lattice)) for _ in 1:9]
        push!(lattices, square_lattice) # add the original lattice as well

        walkers = [LatticeWalker(lattice) for lattice in lattices]
        
        # Create a simple classical hamiltonian for testing
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

        # Test energy assignment without perturbation
        ls = LatticeGasWalkers(deepcopy(walkers), ham; assign_energy=true)

        @test length(ls.walkers) == 10
        @test all(w.energy != 0.0u"eV" for w in ls.walkers)
        @test all(typeof(w.energy) == typeof(1.0u"eV") for w in ls.walkers)
        @test ls.walkers[end].energy == -0.2u"eV" # known value for the original lattice

        # Test energy assignment with perturbation
        ls_perturbed = LatticeGasWalkers(deepcopy(walkers), ham; assign_energy=true, perturb_energy=0.1)
        @test length(ls_perturbed.walkers) == 10
        @test all(w.energy != 0.0u"eV" for w in ls_perturbed.walkers)
        @test all(typeof(w.energy) == typeof(1.0u"eV") for w in ls_perturbed.walkers)
        @test ls_perturbed.walkers[end].energy != ls.walkers[end].energy # should differ due to perturbation
        @test all(abs(ls_perturbed.walkers[i].energy - ls.walkers[i].energy) ≤ 0.1u"eV" for i in 1:10)
    end

end