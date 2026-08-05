# Custom-Hamiltonian fixtures for the extension-contract testset at the
# bottom of this file; type definitions must sit at the file's top level.
# ExtCoordHam: a coordination-dependent global functional no shipped pair,
# cluster, or site-field type expresses (the documentation page's example).
struct ExtCoordHam <: ClassicalHamiltonian
    j::Float64
end
function FreeBird.EnergyEval.interacting_energy(lattice::SLattice, h::ExtCoordHam)
    occ = lattice.components[1]
    doubly = count(s -> occ[s] &&
                        count(n -> occ[n], lattice.neighbors[s][1]) == 2,
                   eachindex(occ))
    return h.j * doubly^2 * u"eV"
end
# ExtDelegatingHam: wraps a shipped Hamiltonian unchanged, so a same-seed
# run must be bit-identical to the shipped type (stream neutrality of the
# extension path, including the setup validation)
struct ExtDelegatingHam <: ClassicalHamiltonian
    inner
end
FreeBird.EnergyEval.interacting_energy(lattice::SLattice, h::ExtDelegatingHam) =
    interacting_energy(lattice, h.inner)
# ExtBadHam: violates the return-type contract (plain Float64)
struct ExtBadHam <: ClassicalHamiltonian end
FreeBird.EnergyEval.interacting_energy(lattice::SLattice, h::ExtBadHam) = 1.6

@testset "AbstractLiveSets Tests" begin
    @testset "atomistic_livesets.jl tests" begin

        # Create periodic box and LJ parameters
        box = [[10.0u"Å", 0u"Å", 0u"Å"],
               [0u"Å", 10.0u"Å", 0u"Å"],
               [0u"Å", 0u"Å", 10.0u"Å"]]
        lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

        ljs = CompositeParameterSets(3, [lj for _ in 1:6]) # 3 components, 6 parameters

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

    @testset "custom-Hamiltonian extension contract" begin
        using Random
        ext_sq = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:16]],
            adsorptions=:full)
        ext_tri = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            supercell_dimensions=(3, 3, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:18]],
            adsorptions=:full)
        ext_ham = ExtCoordHam(0.1)

        # Negative path: the return-type contract violation raises the
        # descriptive ArgumentError at walker setup, naming the method and
        # the required dimension, both at construction and through
        # exact_enumeration (which constructs a liveset internally)
        bad_walkers = [LatticeWalker(deepcopy(ext_sq), energy=0.0u"eV", iter=0)]
        @test_throws ArgumentError LatticeGasWalkers(bad_walkers, ExtBadHam())
        err = try
            LatticeGasWalkers(bad_walkers, ExtBadHam())
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("interacting_energy", err.msg)
        @test occursin("energy-dimensioned", err.msg)
        @test occursin("Float64", err.msg)
        @test_throws ArgumentError exact_enumeration(deepcopy(ext_sq), ExtBadHam())

        # Positive path: one interacting_energy method serves construction
        # on both geometries, with energies matching an independent
        # re-evaluation of the same rule
        Random.seed!(4661)
        function ext_independent(lat)
            occ = lat.components[1]
            c = 0
            for s in eachindex(occ)
                occ[s] || continue
                nb = count(n -> occ[n], lat.neighbors[s][1])
                c += (nb == 2) ? 1 : 0
            end
            return 0.1 * c^2
        end
        for lat in (ext_sq, ext_tri)
            walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                       for _ in 1:6]
            for w in walkers
                w.configuration.components[1] .= rand(Bool, length(w.configuration.components[1]))
            end
            ls = LatticeGasWalkers(walkers, ext_ham)
            for w in ls.walkers
                @test w.energy.val == ext_independent(w.configuration)
            end
        end

        # Positive path, sampled: a short seeded grand-canonical run drives
        # the custom Hamiltonian through the full loop
        Random.seed!(4662)
        walkers = [LatticeWalker(deepcopy(ext_sq), energy=0.0u"eV", iter=0)
                   for _ in 1:8]
        ls = LatticeGasWalkers(walkers, ext_ham; assign_energy=false)
        gc = GrandCanonicalNestedSamplingParameters(mc_steps=20,
            chemical_potential=-0.05, energy_perturbation=1e-9)
        ext_save = SaveEveryN("t_ext.csv", "t_ext.traj", "t_ext.ls",
                              1000000, 1000000, 1000000)
        df, _, _ = grand_canonical_nested_sampling(ls, gc, Int64(100),
            MCGrandCanonicalMoves(), ext_save;
            observables=[:n_occ => (cfg -> Float64(sum(cfg.components[1])))])
        # The 10^6-step save intervals guarantee no t_ext.* file is written
        # in these short runs; the rm calls are defensive cleanup only
        rm.(["t_ext.csv", "t_ext.traj", "t_ext.ls"], force=true)
        @test nrow(df) > 0
        @test all(isfinite, df.energy)
        @test "n_occ" in names(df)
        @test all(isfinite, df.n_occ)
        @test df.n_occ == Float64.(df.num_particles)

        # Positive path, enumerated: the fixed-N spectra match a brute force
        # over all occupations of the independent rule, on the square
        # C(16, 4) and the triangular C(18, 3) sectors
        for (lat0, M, N) in [(ext_sq, 16, 4), (ext_tri, 18, 3)]
            latN = deepcopy(lat0)
            latN.components[1] .= vcat(fill(true, N), fill(false, M - N))
            dfe, _ = exact_enumeration(latN, ext_ham)
            @test nrow(dfe) == binomial(M, N)
            probe = deepcopy(lat0)
            brute = Float64[]
            for mask in 0:(2^M-1)
                count_ones(mask) == N || continue
                for s in 1:M
                    probe.components[1][s] = ((mask >> (s - 1)) & 1) == 1
                end
                push!(brute, ext_independent(probe))
            end
            @test sort([ustrip(u"eV", e) for e in dfe.energy]) == sort(brute)
        end

        # Shipped types still construct through the validated path with
        # bit-identical energies (the validation is a no-op for them), and
        # a delegating custom Hamiltonian is a same-seed drop-in for the
        # shipped type it wraps: identical ideal-gas-referenced ledgers
        ext_gham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        function ext_igref(seed, h)
            Random.seed!(seed)
            ws = [LatticeWalker(deepcopy(ext_sq), energy=0.0u"eV", iter=0)
                  for _ in 1:12]
            l = LatticeGasWalkers(ws, h; assign_energy=false)
            p = IdealGasReferencedGCNSParameters(mc_steps=20,
                reference_fugacity=1.0, energy_perturbation=1e-9)
            d, fl, _ = ideal_gas_referenced_nested_sampling(l, p, Int64(80),
                MCGrandCanonicalMoves(), ext_save)
            rm.(["t_ext.csv", "t_ext.traj", "t_ext.ls"], force=true)
            return d, [w.energy.val for w in fl.walkers]
        end
        dfS, liveS = ext_igref(4663, ext_gham)
        dfD, liveD = ext_igref(4663, ExtDelegatingHam(ext_gham))
        @test dfS.iter == dfD.iter
        @test dfS.emax == dfD.emax
        @test dfS.num_particles == dfD.num_particles
        @test liveS == liveD
    end

end