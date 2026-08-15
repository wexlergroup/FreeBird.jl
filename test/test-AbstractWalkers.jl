@testset "AbstractWalkers Tests" begin
    @testset "atomistic_walkers.jl test" begin
        @testset "AtomWalker struct and functions tests" begin
            
            # Create a test configuration
            at = FreeBirdIO.generate_multi_type_random_starting_config(
                10.0, 
                [2,1,3,4,5,6], 
                particle_types=[:H,:O,:H,:Fe,:Au,:Cl]
            )

            @testset "Constructor tests" begin
                
                # Test basic constructor with default values
                @test_nowarn AtomWalker{6}(at, list_num_par=[2,1,3,4,5,6])
                walker = AtomWalker{6}(at, list_num_par=[2,1,3,4,5,6])

                @test walker.energy == 0.0u"eV"
                @test walker.iter == 0
                @test walker.list_num_par == [2,1,3,4,5,6]
                @test length(walker.frozen) == 6
                @test walker.energy_frozen_part == 0.0u"eV"

                # Test constructor with custom values
                custom_walker = AtomWalker{6}(
                    at,
                    energy = 1.5u"eV",
                    iter = 10,
                    list_num_par = [2,1,3,4,5,6],
                    frozen = [true, false, true, false, false, false],
                    energy_frozen_part = 0.5u"eV"
                )
                
                @test custom_walker.energy == 1.5u"eV"
                @test custom_walker.iter == 10
                @test custom_walker.list_num_par == [2,1,3,4,5,6]
                @test custom_walker.frozen == [true, false, true, false, false, false]
                @test custom_walker.energy_frozen_part == 0.5u"eV"
            end

            @testset "Error handling test" begin

                # Test incorrect list_num_par length
                @test_throws ArgumentError AtomWalker{6}(
                    at,
                    list_num_par=[1,2,3]
                )

                # Test incorrect frozen length
                @test_throws ArgumentError AtomWalker{6}(
                    at,
                    list_num_par=[2,1,3,4,5,6],
                    frozen=[true, false]
                )

                # Test incompatible list_num_par sum
                @test_throws ArgumentError AtomWalker{6}(
                    at,
                    list_num_par=[1,1,1,1,1,1]  # Sum doesn't match configuration length
                )

                # Test default zeros list_num_par (should throw error)
                @test_throws ArgumentError AtomWalker{6}(at)                
            end

            @testset "Convenience constructor tests" begin

                # Test with no frozen species
                @test_nowarn AtomWalker(at)
                walker = AtomWalker(at)
                @test !any(walker.frozen)
                @test sum(walker.list_num_par) == length(at)

                # Test with frozen species
                walker_frozen = AtomWalker(at, freeze_species=[:H])
                @test any(walker_frozen.frozen)
                @test sum(walker_frozen.list_num_par) == length(at)

                # Test merge_same_species
                walker_merged = AtomWalker(at, merge_same_species=true)
                walker_unmerged = AtomWalker(at, merge_same_species=false)
                @test length(walker_merged.list_num_par) < length(walker_unmerged.list_num_par)
                @test sum(walker_merged.list_num_par) == length(at)
                @test sum(walker_unmerged.list_num_par) == length(at)

                # Test with invalid freeze species
                @test_throws ArgumentError AtomWalker(at, freeze_species=[:Si])

                # Test FlexibleSystem
                flex_at = FlexibleSystem(at)
                @test_nowarn AtomWalker(flex_at)
                walker = AtomWalker(flex_at)
                @test !any(walker.frozen)
                @test sum(walker.list_num_par) == length(flex_at)
            end

            @testset "Show method tests" begin
                walker = AtomWalker(at)  # Use the convenience constructor
                
                # Test single walker show
                buf = IOBuffer()
                @test_nowarn show(buf, walker)
                output = String(take!(buf))
                @test contains(output, "AtomWalker")
                @test contains(output, "configuration")
                @test contains(output, "energy")
                @test contains(output, "iter")
                @test contains(output, "list_num_par")
                @test contains(output, "frozen")
                @test contains(output, "energy_frozen_part")
        
                # Test vector of walkers show
                walkers = [walker for _ in 1:3]
                buf = IOBuffer()
                @test_nowarn show(buf, walkers)
                output = String(take!(buf))
                @test contains(output, "Vector{AtomWalker")
                @test contains(output, "[1]")
                @test contains(output, "[2]")
                @test contains(output, "[3]")               
            end
        
            @testset "Single component special case tests" begin
                # Create a single component configuration
                single_config = FreeBirdIO.generate_multi_type_random_starting_config(
                    10.0,
                    [5],
                    particle_types=[:H]
                )
                
                # Test C=1 case
                walker = AtomWalker{1}(single_config)
                @test length(walker.list_num_par) == 1
                @test walker.list_num_par[1] == length(single_config)
                @test length(walker.frozen) == 1
            end
        end


        @testset "update_walker function tests" begin

            # Create a test configuration
            at = FreeBirdIO.generate_multi_type_random_starting_config(
                10.0, 
                [2,1,3,4,5,6], 
                particle_types=[:H,:O,:H,:Fe,:Au,:Cl]
            )

            walker = AtomWalker(at)
        
            @testset "Valid property updates" begin
                # walker = AtomWalker(at)
                initial_config = walker.configuration
        
                # Test updating energy
                @test_nowarn update_walker!(walker, :energy, 1.5u"eV")
                @test walker.energy == 1.5u"eV"
        
                # Test updating iteration count
                @test_nowarn update_walker!(walker, :iter, 10)
                @test walker.iter == 10
        
                # Test updating frozen status
                new_frozen = [true, false, true]
                @test_nowarn update_walker!(walker, :frozen, new_frozen)
                @test walker.frozen == new_frozen
        
                # Test updating energy_frozen_part
                @test_nowarn update_walker!(walker, :energy_frozen_part, 0.5u"eV")
                @test walker.energy_frozen_part == 0.5u"eV"
        
                # Test updating list_num_par
                new_list_num_par = copy(walker.list_num_par)
                @test_nowarn update_walker!(walker, :list_num_par, new_list_num_par)
                @test walker.list_num_par == new_list_num_par
            end
        
            @testset "Type consistency" begin
                # walker = create_test_walker()
        
                # Test energy must have correct units
                @test_throws DimensionError update_walker!(walker, :energy, 1.5)  # No units
                @test_throws DimensionError update_walker!(walker, :energy, 1.5u"m")  # Wrong units
        
                # Test iter must be Integer
                @test_throws InexactError update_walker!(walker, :iter, 10.5)
        
                # Test energy_frozen_part must have correct units
                @test_throws DimensionError update_walker!(walker, :energy_frozen_part, 0.5)  # No units
                @test_throws DimensionError update_walker!(walker, :energy_frozen_part, 0.5u"m")  # Wrong units
        
                # Test list_num_par must be Vector{Int64}
                @test_throws InexactError update_walker!(walker, :list_num_par, [1.5, 2.5, 3.5])
            end
        
            @testset "Invalid property updates" begin
        
                # Test updating non-existent property
                @test_throws ErrorException update_walker!(walker, :nonexistent_property, 42)
                
                # Test with invalid key type
                @test_throws MethodError update_walker!(walker, "energy", 1.5u"eV")  # String instead of Symbol
            end
        
            @testset "Return value" begin
                
                # Test that the function returns the updated walker
                result = update_walker!(walker, :energy, 2.0u"eV")
                @test result === walker  # Test that it returns the same object (not a copy)
                @test result.energy == 2.0u"eV"  # Test that the update was applied
            end
        
            @testset "Multiple updates" begin
                
                # Test multiple sequential updates
                @test_nowarn begin
                    update_walker!(walker, :energy, 1.0u"eV")
                    update_walker!(walker, :iter, 5)
                    update_walker!(walker, :energy_frozen_part, 0.3u"eV")
                end
                
                # Verify all updates were applied
                @test walker.energy == 1.0u"eV"
                @test walker.iter == 5
                @test walker.energy_frozen_part == 0.3u"eV"
            end
        end
    end


    @testset "lattice_walkers.jl tests" begin
        @testset "compute_neighbors function tests" begin

            # Create a simple cubic lattice
            lattice_vectors = [
                2.0 0.0 0.0;
                0.0 2.0 0.0;
                0.0 0.0 2.0
            ]

            @testset "Non-periodic tests" begin
                positions = [
                    0.0 0.0 0.0;  # atom 1 at origin
                    1.0 0.0 0.0;  # atom 2
                    0.0 1.0 0.0;  # atom 3
                    1.0 1.0 0.0;  # atom 4
                    0.0 0.0 1.0;  # atom 5
                    1.0 0.0 1.0;  # atom 6
                    0.0 1.0 1.0;  # atom 7
                    1.0 1.0 1.0   # atom 8
                    ]

                periodicity = (false, false, false)
                cutoff_radii = [1.1, 1.8]  # First and second nearest neighbors
                
                neighbors = AbstractWalkers.compute_neighbors(lattice_vectors, positions, periodicity, cutoff_radii)
                
                # Check atom 1 (corner atom)
                @test sort(neighbors[1][1]) == [2, 3, 5]
                @test sort(neighbors[1][2]) == [4, 6, 7, 8]
                
                # Check atom 8 (opposite corner)
                @test sort(neighbors[8][1]) == [4, 6, 7]  # First neighbors
                @test sort(neighbors[8][2]) == [1, 2, 3, 5]  # Second neighbors
            end
                    
            @testset "Fully periodic tests" begin

                positions = [
                    0.0 0.0 0.0;  # atom 1 at origin
                    1.0 1.0 0.0;  # atom 2
                    1.0 0.0 1.0;  # atom 3
                    0.0 1.0 1.0;  # atom 4
                    1.0 1.0 1.0   # atom 5
                    ]

                periodicity = (true, true, true)
                cutoff_radii = [1.1, 1.5, 1.8]

                # Every cutoff exceeds half the cell edge on this one-cell
                # torus, so pairs are connected through several images: the
                # pinned per-site counts below are the documented
                # minimum-image convention, and the collapsed-image warning
                # fires
                neighbors = @test_logs (:warn, r"wrap the periodic cell") match_mode=:any AbstractWalkers.compute_neighbors(
                    lattice_vectors, positions, periodicity, cutoff_radii)
                
                for i in 1:5
                    if i == 1
                        @test length(neighbors[i][1]) == 0
                        @test length(neighbors[i][2]) == 3
                        @test length(neighbors[i][3]) == 1
                    elseif i == 5
                        @test length(neighbors[i][1]) == 3
                        @test length(neighbors[i][2]) == 0
                        @test length(neighbors[i][3]) == 1
                    else
                        @test length(neighbors[i][1]) == 1
                        @test length(neighbors[i][2]) == 3
                        @test length(neighbors[i][3]) == 0
                    end
                end
            end
        end


        @testset "image multiplicity and cutoff_radii validation" begin

            # Full-occupancy single-component square builder
            sq_lattice(dims, per, cutoffs; kwargs...) = MLattice{1,SquareLattice}(;
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=dims,
                periodicity=per,
                cutoff_radii=cutoffs,
                components=[[true for _ in 1:prod(dims)]],
                adsorptions=:full,
                kwargs...)

            @testset "wrapped 4x4 three-shell coordination pins" begin
                # The third shell (distance 2 = L/2) wraps: the +2 and -2
                # images of the same site collapse to one entry under the
                # minimum-image convention ...
                lat_def = @test_logs (:warn, r"wrap the periodic cell") match_mode=:any sq_lattice(
                    (4, 4, 1), (true, true, false), [1.1, 1.5, 2.1])
                @test all(length.(lat_def.neighbors[i]) == [4, 4, 2] for i in 1:16)
                # ... while multiplicity restores the bulk-tiled coordination,
                # with no warning
                lat_mult = @test_logs min_level=Base.CoreLogging.Warn sq_lattice(
                    (4, 4, 1), (true, true, false), [1.1, 1.5, 2.1], image_multiplicity=true)
                @test all(length.(lat_mult.neighbors[i]) == [4, 4, 4] for i in 1:16)
                # Element-pinned lists for the corner site: entries ascending
                # in j, multiplicity entries grouped by j
                @test lat_def.neighbors[1] == [[2, 4, 5, 13], [6, 8, 14, 16], [3, 9]]
                @test lat_mult.neighbors[1] == [[2, 4, 5, 13], [6, 8, 14, 16], [3, 3, 9, 9]]
            end

            @testset "L = 2r degeneracy and warning content: triangular (4, 2, 1)" begin
                # The a2 circumference 2√3 equals exactly twice the
                # second-neighbor distance √3, so the degenerate image
                # distances are bit-identical — no tolerance involved — and
                # the warning names the wrapped shell with its exact
                # collapsed-image count
                warnpat = r"neighbor shell\(s\) \[2 \(cutoff 1\.8\): 16 collapsed image bond\(s\)\] wrap the periodic cell"
                lat_def = @test_logs (:warn, warnpat) match_mode=:any MLattice{1,TriangularLattice}()
                M = length(lat_def.components[1])
                @test M == 16
                @test all(length(lat_def.neighbors[i][1]) == 6 for i in 1:M)
                @test all(length(lat_def.neighbors[i][2]) == 5 for i in 1:M)

                lat_mult = @test_logs min_level=Base.CoreLogging.Warn MLattice{1,TriangularLattice}(
                    image_multiplicity=true)
                @test all(length(lat_mult.neighbors[i][1]) == 6 for i in 1:M)
                @test all(length(lat_mult.neighbors[i][2]) == 6 for i in 1:M)
                # Exactly one second-shell neighbor per site — the a2-wrapped
                # one — is doubled
                @test all(length(unique(lat_mult.neighbors[i][2])) == 5 for i in 1:M)
                # The faithful first shell agrees between the conventions
                @test all(lat_mult.neighbors[i][1] == lat_def.neighbors[i][1] for i in 1:M)
            end

            @testset "1D chains, periodicity (true, false, false)" begin
                # (4,1,1): the second shell (distance 2 = L/2) wraps
                ch_def = @test_logs (:warn, r"wrap the periodic cell") match_mode=:any sq_lattice(
                    (4, 1, 1), (true, false, false), [1.1, 2.1])
                @test [n[1] for n in ch_def.neighbors] == [[2, 4], [1, 3], [2, 4], [1, 3]]
                @test [n[2] for n in ch_def.neighbors] == [[3], [4], [1], [2]]
                ch_mult = @test_logs min_level=Base.CoreLogging.Warn sq_lattice(
                    (4, 1, 1), (true, false, false), [1.1, 2.1], image_multiplicity=true)
                @test [n[1] for n in ch_mult.neighbors] == [[2, 4], [1, 3], [2, 4], [1, 3]]
                @test [n[2] for n in ch_mult.neighbors] == [[3, 3], [4, 4], [1, 1], [2, 2]]

                # (2,1,1): first-shell multiplicity 2 through the ±1 images of
                # the other site, plus two second-shell self-image entries
                # (j == i) through the ±2 self-images
                ch2_def = @test_logs (:warn, r"wrap the periodic cell") match_mode=:any sq_lattice(
                    (2, 1, 1), (true, false, false), [1.1, 2.1])
                @test ch2_def.neighbors == [[[2], Int[]], [[1], Int[]]]
                ch2_mult = @test_logs min_level=Base.CoreLogging.Warn sq_lattice(
                    (2, 1, 1), (true, false, false), [1.1, 2.1], image_multiplicity=true)
                @test ch2_mult.neighbors == [[[2, 2], [1, 1]], [[1, 1], [2, 2]]]
            end

            @testset "ArgumentError on invalid cutoff_radii ladders" begin
                lv = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 1.0]
                pos = [0.0 0.0 0.0; 1.0 0.0 0.0]
                per = (true, true, false)
                for bad in ([1.5, 1.1], [1.1, 1.1], [-1.0], Float64[], [Inf], [1.1, Inf], [NaN])
                    @test_throws ArgumentError AbstractWalkers.compute_neighbors(lv, pos, per, bad)
                end
                err = try
                    AbstractWalkers.compute_neighbors(lv, pos, per, [1.5, 1.1])
                catch e
                    e
                end
                @test occursin("nested cutoff ladder", err.msg)
                # ... and the validation is reached through the constructors
                @test_throws ArgumentError MLattice{1,SquareLattice}(cutoff_radii=[1.5, 1.1])
                @test_throws ArgumentError MLattice{1,TriangularLattice}(cutoff_radii=Float64[])
            end
        end


        @testset "lattice_positions function tests" begin  
            @testset "Simple cubic lattice" begin

                # Create a simple cubic lattice
                lattice_vectors = [
                    2.0 0.0 0.0;
                    0.0 2.0 0.0;
                    0.0 0.0 2.0
                ]

                basis = [(0.0, 0.0, 0.0)]

                @testset "1x1x1 supercell" begin
                    dimensions = (1, 1, 1)
                    positions = AbstractWalkers.lattice_positions(lattice_vectors, basis, dimensions)
                    
                    @test size(positions) == (1, 3)
                    @test positions[1, :] ≈ [0.0, 0.0, 0.0]
                end

                @testset "2x2x2 supercell" begin
                    dimensions = (2, 2, 2)
                    positions = AbstractWalkers.lattice_positions(lattice_vectors, basis, dimensions)
                    
                    @test size(positions) == (8, 3)
                    expected_positions = [
                        0.0 0.0 0.0;
                        2.0 0.0 0.0;
                        0.0 2.0 0.0;
                        2.0 2.0 0.0;
                        0.0 0.0 2.0;
                        2.0 0.0 2.0;
                        0.0 2.0 2.0;
                        2.0 2.0 2.0
                    ]
                    @test positions ≈ expected_positions
                end
            end

            @testset "Body-centered cubic" begin

                lattice_vectors = [
                    2.0 0.0 0.0;
                    0.0 2.0 0.0;
                    0.0 0.0 2.0
                ]

                basis = [
                    (0.0, 0.0, 0.0),
                    (1.0, 1.0, 1.0)
                ]
                
                # Test 1x1x1 cell
                dimensions = (1, 1, 1)
                positions = AbstractWalkers.lattice_positions(lattice_vectors, basis, dimensions)
                
                @test size(positions) == (2, 3)
                @test positions[1, :] ≈ [0.0, 0.0, 0.0]
                @test positions[2, :] ≈ [1.0, 1.0, 1.0]
            end

            @testset "Face-centered cubic" begin

                lattice_vectors = [
                    2.0 0.0 0.0;
                    0.0 2.0 0.0;
                    0.0 0.0 2.0
                ]

                basis = [
                    (0.0, 0.0, 0.0),
                    (1.0, 1.0, 0.0),
                    (1.0, 0.0, 1.0),
                    (0.0, 1.0, 1.0)
                ]
                
                # Test 1x1x1 cell
                dimensions = (1, 1, 1)
                positions = AbstractWalkers.lattice_positions(lattice_vectors, basis, dimensions)
                
                @test size(positions) == (4, 3)
                @test positions[1, :] ≈ [0.0, 0.0, 0.0]
                @test positions[2, :] ≈ [1.0, 1.0, 0.0]
                @test positions[3, :] ≈ [1.0, 0.0, 1.0]
                @test positions[4, :] ≈ [0.0, 1.0, 1.0]
            end

            @testset "Complex basis" begin
                lattice_vectors = [
                    2.0 0.0 0.0;
                    0.0 2.0 0.0;
                    0.0 0.0 2.0
                ]

                basis = [
                    (0.0, 0.0, 0.0),
                    (0.25, 0.25, 0.25),
                    (0.75, 0.75, 0.75)
                ]
                
                dimensions = (1, 1, 1)
                positions = AbstractWalkers.lattice_positions(lattice_vectors, basis, dimensions)
                
                @test size(positions) == (3, 3)
                @test positions[1, :] ≈ [0.0, 0.0, 0.0]
                @test positions[2, :] ≈ [0.25, 0.25, 0.25]
                @test positions[3, :] ≈ [0.75, 0.75, 0.75]
            end
        end


        @testset "split_into_subarrays function tests" begin

            @testset "Error Handling" begin
                # Test with N ≤ 0
                @test_throws DivideError AbstractWalkers.split_into_subarrays([1, 2, 3], 0)
            end

            @testset "Even Division" begin
                # Test array with length divisible by N
                arr = collect(1:9)
                result = AbstractWalkers.split_into_subarrays(arr, 3)
                @test length(result) == 3
                @test result == [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
            end

            @testset "Uneven Division" begin
                # Test array with length not divisible by N
                arr = collect(1:10)
                result = AbstractWalkers.split_into_subarrays(arr, 3)
                @test length(result) == 3
                @test result == [[1, 2, 3, 4], [5, 6, 7], [8, 9, 10]]
                
                # Test when remainder is distributed
                arr = collect(1:7)
                result = AbstractWalkers.split_into_subarrays(arr, 3)
                @test length(result) == 3
                @test result == [[1, 2, 3], [4, 5], [6, 7]]
                
                # Verify the first subarrays get the extra elements
                arr = collect(1:8)
                result = AbstractWalkers.split_into_subarrays(arr, 3)
                @test length(result[1]) == 3  # First subarray gets extra element
                @test length(result[2]) == 3  # Second subarray gets extra element
                @test length(result[3]) == 2  # Last subarray gets base size
            end

            @testset "Edge Cases" begin
                # Test empty array
                result = AbstractWalkers.split_into_subarrays(Int[], 3)
                @test length(result) == 3
                @test all(isempty, result)
                
                # Test single element array
                arr = [1]
                result = AbstractWalkers.split_into_subarrays(arr, 3)
                @test length(result) == 3
                @test result == [[1], Int[], Int[]]
                
                # Test when N equals array length
                arr = collect(1:4)
                result = AbstractWalkers.split_into_subarrays(arr, 4)
                @test length(result) == 4
                @test all(x -> length(x) == 1, result)
                @test result == [[1], [2], [3], [4]]
                
                # Test when N is larger than array length
                arr = collect(1:3)
                result = AbstractWalkers.split_into_subarrays(arr, 5)
                @test length(result) == 5
                @test count(!isempty, result) == 3
                @test result == [[1], [2], [3], Int[], Int[]]
            end

            @testset "Different Types" begin
                # Test with Float64
                arr = [1.5, 2.5, 3.5, 4.5]
                result = AbstractWalkers.split_into_subarrays(arr, 2)
                @test eltype(result[1]) == Float64
                @test result == [[1.5, 2.5], [3.5, 4.5]]
                
                # Test with String
                arr = ["a", "b", "c", "d", "e"]
                result = AbstractWalkers.split_into_subarrays(arr, 3)
                @test eltype(result[1]) == String
                @test result == [["a", "b"], ["c", "d"], ["e"]]
                
                # Test with custom type
                struct CustomType
                    value::Int
                end
                arr = [CustomType(i) for i in 1:4]
                result = AbstractWalkers.split_into_subarrays(arr, 2)
                @test eltype(result[1]) == CustomType
                @test length(result) == 2
                @test length(result[1]) == 2
                @test length(result[2]) == 2
            end
        end


        @testset "mlattice_setup function tests" begin
            
            # Common test parameters
            basis_single = [(0.0, 0.0, 0.0)]
            basis_double = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
            dims_2d = (2, 2, 1)
            
            components_cum1 = [[true, false, false, false], [false, false, false, true]]
            adsorptions_cum1 = [false, true, true, false]
            components_cum2 = [[1], [4]]
            adsorptions_cum2 = [2, 3]
            
            @testset "Equal components and Full adsorption" begin
                # Test with single basis point
                components, adsorptions = AbstractWalkers.mlattice_setup(
                    2,
                    basis_single,
                    dims_2d,
                    :equal,
                    :full
                )
                
                @test length(components) == 2  # Two components
                @test length(components[1]) == 4  # 2×2×1 = 4 sites
                @test count(components[1]) == 2  # Half of sites in first component
                @test count(components[2]) == 2  # Half of sites in second component
                @test all(adsorptions)  # All sites are adsorption sites
                
                # Test with double basis points
                components, adsorptions = AbstractWalkers.mlattice_setup(
                    3,  # 3 components
                    basis_double,
                    dims_2d,
                    :equal,
                    :full
                )
                
                @test length(components) == 3  # Three components
                @test length(components[1]) == 8  # 2×2×1×2 = 8 sites
                @test all(count.(components) .≈ [3, 3, 2])  # Distribution of sites

                @testset "Custom components and adsorptions" begin

                    # components isa Vector{Vector{Bool}} and adsoprtions isa Vector{Bool}
                    components, adsorptions = AbstractWalkers.mlattice_setup(
                        2,
                        basis_single,
                        dims_2d,
                        components_cum1,
                        adsorptions_cum1
                    )

                    @test length(components) == 2  
                    @test length(components[1]) == 4  
                    @test count(components[1]) == 1  
                    @test count(components[2]) == 1

                    # components isa Vector{Vector{Int}} and adsoprtions isa Vector{Int}
                    components, adsorptions = AbstractWalkers.mlattice_setup(
                        2,
                        basis_single,
                        dims_2d,
                        components_cum2,
                        adsorptions_cum2
                    )

                    @test length(components) == 2  
                    @test length(components[1]) == 4  
                    @test count(components[1]) == 1  
                    @test count(components[2]) == 1

                    # Error cases for components and adsoprtions
                    @test_throws ArgumentError AbstractWalkers.mlattice_setup(
                        2,
                        basis_single,
                        dims_2d,
                        :not_equal,
                        :not_full
                    )
                end
            end
        end


        @testset "MLattice struct and funcitons tests" begin
            @testset "Square lattice constructor" begin

                # Default parameters
                lattice = MLattice{2,SquareLattice}()

                @test size(lattice.lattice_vectors) == (3, 3)
                @test lattice.lattice_vectors ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
                @test lattice.basis == [(0.0, 0.0, 0.0)]
                @test lattice.supercell_dimensions == (4, 4, 1)
                @test lattice.periodicity == (true, true, false)
                @test lattice.cutoff_radii == [1.1, 1.5]
                @test length(lattice.neighbors) == 16
                @test length(lattice.components) == 2
                @test length(lattice.adsorptions) == 16
                @test all(lattice.adsorptions)

                # Custom parameters
                custom_basis = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
                custom_dims = (2, 2, 1)
                custom_periodicity = (false, false, false)
                custom_components = [[1, 2], [3, 4]]
                custom_adsorptions = [1, 3]

                lattice = MLattice{2,SquareLattice}(
                    lattice_constant=2.0,
                    basis=custom_basis,
                    supercell_dimensions=custom_dims,
                    periodicity=custom_periodicity,
                    components=custom_components,
                    adsorptions=custom_adsorptions
                )
                
                # Isotropic default: the stored third diagonal now follows
                # lattice_constant (at d3 = 1 the z axis never enters
                # positions, so the geometry is unchanged)
                @test lattice.lattice_vectors ≈ [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 2.0]
                @test lattice.basis == custom_basis
                @test lattice.supercell_dimensions == custom_dims
                @test lattice.periodicity == custom_periodicity
                @test length(lattice.neighbors) == 8
                @test length(lattice.components) == length(custom_components)
                @test count(lattice.adsorptions) == 2
            end

            @testset "Triangular lattice constructor" begin

                # Default parameters
                lattice = MLattice{2,TriangularLattice}()
            
                @test size(lattice.lattice_vectors) == (3, 3)
                @test lattice.lattice_vectors ≈ [1.0 0.0 0.0; 0.0 sqrt(3) 0.0; 0.0 0.0 1.0]
                @test lattice.basis == [(0.0, 0.0, 0.0), (1/2, sqrt(3)/2, 0.0)]
                @test lattice.supercell_dimensions == (4, 2, 1)
                @test lattice.periodicity == (true, true, false)
                # Triangular default is [1.1, 1.8]: the second-neighbor
                # distance √3 ≈ 1.732 lies beyond the square-lattice 1.5
                @test lattice.cutoff_radii == [1.1, 1.8]
                @test length(lattice.neighbors) == 16
                @test length(lattice.components) == 2
                @test length(lattice.adsorptions) == 16
                @test all(lattice.adsorptions)

                # Custom parameters
                custom_basis = [(0.0, 0.0, 0.0)]
                custom_dims = (2, 2, 1)
                custom_periodicity = (false, false, false)
                custom_components = [[1, 2], [3, 4]]
                custom_adsorptions = [1, 3]
                
                lattice = MLattice{2,TriangularLattice}(
                    lattice_constant=2.0,
                    basis=custom_basis,
                    supercell_dimensions=custom_dims,
                    periodicity=custom_periodicity,
                    components=custom_components,
                    adsorptions=custom_adsorptions
                )
            
                # Isotropic default: see the square-constructor twin above
                @test lattice.lattice_vectors ≈ [2.0 0.0 0.0; 0.0 2.0*sqrt(3) 0.0; 0.0 0.0 2.0]
                @test lattice.basis == custom_basis
                @test lattice.supercell_dimensions == custom_dims
                @test lattice.periodicity == custom_periodicity
                @test length(lattice.neighbors) == 4
                @test length(lattice.components) == 2
                @test count(lattice.adsorptions) == 2
            end

            @testset "Interlayer spacing" begin
                # Isotropic fix: with lattice_constant = 2.0 the out-of-plane
                # spacing now follows it, so the nearest-neighbor shell mixes
                # in-plane and interlayer bonds at the same distance 2.0
                # (previously the interlayer bonds sat at 1.0, silently)
                iso = MLattice{1,SquareLattice}(
                    lattice_constant=2.0,
                    supercell_dimensions=(3, 3, 2),
                    cutoff_radii=[2.2],
                    components=[[false for _ in 1:18]]
                )
                @test iso.lattice_vectors[3, 3] == 2.0
                # lattice_positions ordering: dimension 3 outermost, so the
                # second layer is the second block of 9 sites
                @test iso.positions[:, 3] == vcat(zeros(9), fill(2.0, 9))
                # 4 in-plane + 1 axial on the two-layer slab (z non-periodic)
                @test all(length(iso.neighbors[s][1]) == 5 for s in 1:18)
                @test all((s + 9) in iso.neighbors[s][1] for s in 1:9)

                # Old-geometry reproduction: interlayer_spacing = 1.0 restores
                # the pre-change mixed-scale positions exactly
                old = MLattice{1,SquareLattice}(
                    lattice_constant=2.0,
                    interlayer_spacing=1.0,
                    supercell_dimensions=(3, 3, 2),
                    cutoff_radii=[2.2],
                    components=[[false for _ in 1:18]]
                )
                @test old.lattice_vectors[3, 3] == 1.0
                @test old.positions[:, 3] == vcat(zeros(9), ones(9))

                # Tetragonal shell separation: a = 1.0, c = 1.25 with the
                # ladder [1.1, 1.35] puts in-plane bonds alone in shell 1 and
                # interlayer bonds alone in shell 2 (next distance sqrt(2) ≈
                # 1.414 excluded); construction is silent (in-plane
                # circumference 3 > 2*1.1, no collapsed images; no empty shell)
                tet = @test_logs min_level = Base.CoreLogging.Warn MLattice{1,SquareLattice}(
                    lattice_constant=1.0,
                    interlayer_spacing=1.25,
                    supercell_dimensions=(3, 3, 2),
                    cutoff_radii=[1.1, 1.35],
                    components=[[false for _ in 1:18]]
                )
                @test all(length(tet.neighbors[s][1]) == 4 for s in 1:18)
                @test all(length(tet.neighbors[s][2]) == 1 for s in 1:18)
                @test all(tet.neighbors[s][2] == [s + 9] for s in 1:9)
                # Three layers: the middle layer has axial neighbors both ways
                tet3 = MLattice{1,SquareLattice}(
                    lattice_constant=1.0,
                    interlayer_spacing=1.25,
                    supercell_dimensions=(3, 3, 3),
                    cutoff_radii=[1.1, 1.35],
                    components=[[false for _ in 1:27]]
                )
                @test all(length(tet3.neighbors[s][2]) == 1 for s in 1:9)
                @test all(length(tet3.neighbors[s][2]) == 2 for s in 10:18)
                @test all(length(tet3.neighbors[s][2]) == 1 for s in 19:27)

                # Default-path no-change: an explicit spacing equal to
                # lattice_constant reproduces the default element for element
                expl = MLattice{1,SquareLattice}(
                    lattice_constant=2.0,
                    interlayer_spacing=2.0,
                    supercell_dimensions=(3, 3, 2),
                    cutoff_radii=[2.2],
                    components=[[false for _ in 1:18]]
                )
                @test expl.positions == iso.positions
                @test expl.neighbors == iso.neighbors

                # Guard paths: finite and strictly positive
                @test_throws ArgumentError MLattice{1,SquareLattice}(
                    interlayer_spacing=0.0, components=[[false for _ in 1:16]])
                @test_throws ArgumentError MLattice{1,SquareLattice}(
                    interlayer_spacing=-1.0, components=[[false for _ in 1:16]])
                @test_throws ArgumentError MLattice{1,SquareLattice}(
                    interlayer_spacing=Inf, components=[[false for _ in 1:16]])
                @test_throws ArgumentError MLattice{1,TriangularLattice}(
                    interlayer_spacing=0.0, components=[[false for _ in 1:16]])

                # Triangular acceptance: keyword lands in the third vector
                tri = MLattice{1,TriangularLattice}(
                    supercell_dimensions=(4, 2, 2),
                    interlayer_spacing=1.25,
                    cutoff_radii=[1.1, 1.3],
                    components=[[false for _ in 1:32]]
                )
                @test tri.lattice_vectors ≈ [1.0 0.0 0.0; 0.0 sqrt(3) 0.0; 0.0 0.0 1.25]
            end

            @testset "Error handling" begin
                # Test wrong number of components
                @test_throws ArgumentError MLattice{2,SquareLattice}(
                    lattice_constant=1.0,
                    supercell_dimensions=(2, 2, 1),
                    components=[[true, false, true, false]]
                )

                # Test invalid number of components
                @test_throws BoundsError MLattice{3,SquareLattice}(
                    components=[[1], [2]]  # Only 2 components for 3-component system
                )
            end

            @testset "Component distribution" begin
                # Test equal distribution
                lattice = MLattice{3,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=:equal
                )
                
                component_sums = [sum(comp) for comp in lattice.components]
                @test component_sums ≈ [2, 1, 1]  # 4 sites distributed as evenly as possible
                
                # Test custom distribution with boolean vectors
                bool_components = [[true, false, true, false],
                                 [false, true, false, true]]
                lattice = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=bool_components
                )
                
                @test lattice.components == bool_components
            end

            @testset "Adsorption sites" begin
                # Test full adsorption
                lattice = MLattice{2,SquareLattice}(
                    lattice_constant=1.0,
                    supercell_dimensions=(2, 2, 1),
                    adsorptions=:full
                )
                @test all(lattice.adsorptions)
                
                # Test custom adsorption sites
                custom_adsorptions = [1, 3]  # Indices of adsorption sites
                lattice_custom = MLattice{2,SquareLattice}(
                    lattice_constant=1.0,
                    supercell_dimensions=(2, 2, 1),
                    adsorptions=custom_adsorptions
                )
                # Check if the correct sites are marked as adsorption sites
                expected_adsorptions = [true, false, true, false]
                @test lattice_custom.adsorptions == expected_adsorptions
            end
        end


        @testset "Lattice type aliases and utility functions" begin

            @testset "SLattice type alias" begin
                # Test SLattice with SquareLattice
                square_lattice = SLattice{SquareLattice}()
                @test AbstractWalkers.num_lattice_components(square_lattice) == 1
                @test square_lattice isa MLattice{1,SquareLattice}
                
                # Test SLattice with TriangularLattice
                triangular_lattice = SLattice{TriangularLattice}()
                @test AbstractWalkers.num_lattice_components(triangular_lattice) == 1
                @test triangular_lattice isa MLattice{1,TriangularLattice}
                
                # Test constructor parameters work with alias
                custom_square = SLattice{SquareLattice}(
                    lattice_constant=2.0,
                    supercell_dimensions=(2, 2, 1)
                )

                @test custom_square.lattice_vectors[1,1] ≈ 2.0
                @test custom_square.supercell_dimensions == (2, 2, 1)
            end

            @testset "GLattice type alias" begin
                # Test GLattice
                generic_2comp = GLattice{2}
                @test generic_2comp == MLattice{2,GenericLattice}
            end

            @testset "num_lattice_components function" begin
                # Test with different component numbers
                @test AbstractWalkers.num_lattice_components(MLattice{1,SquareLattice}()) == 1
                @test AbstractWalkers.num_lattice_components(MLattice{2,TriangularLattice}()) == 2
                
                # Test with type aliases
                @test AbstractWalkers.num_lattice_components(SLattice{SquareLattice}()) == 1
            end

            @testset "num_sites function" begin
                # Test square lattice with single basis point
                square_lattice = MLattice{1,SquareLattice}(
                    supercell_dimensions=(2, 3, 1),
                    basis=[(0.0, 0.0, 0.0)]
                )
                @test num_sites(square_lattice) == 6  # 2×3×1×1
        
                # Test triangular lattice with two basis points
                triangular_lattice = MLattice{1,TriangularLattice}(
                    supercell_dimensions=(2, 2, 1),
                    basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
                )
                @test num_sites(triangular_lattice) == 8  # 2×2×1×2
        
                # Test generic lattice with multiple basis points
                generic_lattice = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 2, 2),
                    basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (0.5, 0.0, 0.5)]
                )
                @test num_sites(generic_lattice) == 24  # 2×2×2×3
        
                # Edge cases
                edge_lattice = MLattice{1,SquareLattice}(
                    supercell_dimensions=(1, 1, 1),
                    basis=[(0.0, 0.0, 0.0)]
                )
                @test num_sites(edge_lattice) == 1  # 1×1×1×1
        
                large_lattice = MLattice{1,SquareLattice}(
                    supercell_dimensions=(10, 10, 1),
                    basis=[(0.0, 0.0, 0.0)]
                )
                @test num_sites(large_lattice) == 100  # 10×10×1×1
            end

            @testset "occupied_site_count function" begin
                # Test equal distribution
                lattice_equal = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=:equal
                )
                occupancy_equal = occupied_site_count(lattice_equal)
                @test length(occupancy_equal) == 2
                @test sum(occupancy_equal) == 4  # Total sites
                @test occupancy_equal ≈ [2, 2]  # Equal distribution
                
                # Test custom distribution with components specified by indices
                custom_components = [[1, 2], [3, 4]]
                lattice_custom = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=custom_components
                )
                occupancy_custom = occupied_site_count(lattice_custom)
                @test occupancy_custom == [2, 2]
                
                # Test uneven distribution
                uneven_lattice = MLattice{3,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=:equal
                )
                occupancy_uneven = occupied_site_count(uneven_lattice)
                @test length(occupancy_uneven) == 3
                @test sum(occupancy_uneven) == 4  # Total sites
                @test occupancy_uneven[1] ≥ occupancy_uneven[2] ≥ occupancy_uneven[3]  # Distribution pattern
                
                # Test with empty components
                empty_components = [[true, false, false, false], [false, false, false, false]]
                lattice_empty = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=empty_components
                )
                occupancy_empty = occupied_site_count(lattice_empty)
                @test occupancy_empty == [1, 0]
                
                # Test with fully occupied components
                full_components = [[true, true, true, true], [false, false, false, false]]
                lattice_full = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=full_components
                )
                occupancy_full = occupied_site_count(lattice_full)
                @test occupancy_full == [4, 0]
            end
        end


        @testset "LatticeWalker struct and functions tests" begin
            # Create test lattices
            square_lattice = MLattice{2,SquareLattice}()
            triangular_lattice = MLattice{3,TriangularLattice}()
        
            @testset "Constructor" begin
                # Test default constructor
                walker = LatticeWalker(square_lattice)
                @test walker isa LatticeWalker{2}
                @test walker.configuration === square_lattice
                @test walker.energy == 0.0u"eV"
                @test walker.iter == 0
        
                # Test custom parameters
                custom_walker = LatticeWalker(triangular_lattice, energy=1.5u"eV", iter=10)
                @test custom_walker isa LatticeWalker{3}
                @test custom_walker.energy == 1.5u"eV"
                @test custom_walker.iter == 10
            end
        
            @testset "show methods" begin
                # Test single walker display
                walker = LatticeWalker(square_lattice, energy=2.0u"eV", iter=5)
                output = sprint(show, walker)
                @test contains(output, "LatticeWalker(")
                @test contains(output, "configuration")
                @test contains(output, "energy: 2.0 eV")
                @test contains(output, "iter: 5")
        
                # Test vector of walkers display
                walkers = [
                    LatticeWalker(square_lattice, energy=1.0u"eV", iter=1),
                    LatticeWalker(square_lattice, energy=2.0u"eV", iter=2)
                ]
                output_vec = sprint(show, walkers)
                @test contains(output_vec, "LatticeWalker{2}")
            end
        
            @testset "Mutable behavior" begin
                walker = LatticeWalker(square_lattice)
                walker.energy = 3.5u"eV"
                walker.iter = 15
                @test walker.energy == 3.5u"eV"
                @test walker.iter == 15
            end
        end
    end


    @testset "helpers.jl tests" begin

        @testset "custom_sort function tests" begin

        end

        @testset "MLattice show function tests" begin
            
            lattice = MLattice{2,SquareLattice}()
            output = sprint(show, lattice)

            lattice2 = MLattice{2,SquareLattice}(adsorptions=[1])
            output2 = sprint(show, lattice2)

            # Test presence of all required fields
            @test contains(output, "MLattice{2, SquareLattice}")
            @test contains(output, "lattice_vectors      : [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]")
            @test contains(output, "positions            : 16 grid points")
            @test contains(output, "supercell_dimensions : (4, 4, 1)")
            @test contains(output, "basis                : [(0.0, 0.0, 0.0)]")
            @test contains(output, "periodicity          : (true, true, false)")
            @test contains(output, "cutoff radii         : 2 nearest neighbors")
            @test contains(output, "adsorptions          : full adsorption")

            @test contains(output2, "adsorptions          : \n      ● ○ ○ ○ \n      ○ ○ ○ ○ \n      ○ ○ ○ ○ \n      ○ ○ ○ ○ \n")
        end

        @testset "merge_components function tests" begin
            # Test 2-component square lattice
            lattice2 = MLattice{2,SquareLattice}(
                supercell_dimensions = (2, 2, 1),
                components = [[true, false, true, false], [false, true, false, true]]
            )
            merged = AbstractWalkers.merge_components(lattice2)
            @test merged == [1, 2, 1, 2]
            @test length(merged) == prod(lattice2.supercell_dimensions) * length(lattice2.basis)
    
            # Test 3-component lattice with overlapping components
            lattice3 = MLattice{3,SquareLattice}(
                supercell_dimensions = (2, 2, 1),
                components = [
                    [true, false, false, false],
                    [false, true, false, false],
                    [false, false, true, true]
                ]
            )
            merged3 = AbstractWalkers.merge_components(lattice3)
            @test merged3 == [1, 2, 3, 3]
        end

        @testset "print functions tests" begin
            # Define the general print function
            function capture_print(f, args...)
                io = IOBuffer()
                f(io, args...)
                String(take!(io))
            end

            # Single component print
            s_lattice1 = MLattice{1,SquareLattice}(
                supercell_dimensions = (2, 2, 1),
                components = [[true, false, true, false]]
            )

            s_lattice2 = MLattice{1,SquareLattice}(
                basis= [(1.0, 1.0, 11.0)],
                supercell_dimensions = (2, 2, 2),
                components = [fill(true, 8)]
            )

            # Test print_layer_single_comp()
            print_layer_output1 = capture_print(AbstractWalkers.print_layer, s_lattice1, s_lattice1.components[1])
            @test contains(print_layer_output1, "●")
            @test contains(print_layer_output1, "○")
            @test count("●", print_layer_output1) == 2
            @test count("○", print_layer_output1) == 2

            # Test occupation printing
            occ_output1 = capture_print(AbstractWalkers.print_occupation, s_lattice1)
            occ_output2 = capture_print(AbstractWalkers.print_occupation, s_lattice2)
            @test contains(occ_output1, "●")
            @test contains(occ_output1, "○")
            @test contains(occ_output2, "●")

            # Test adsorption printing
            ads_output1 = capture_print(AbstractWalkers.print_adsorption, s_lattice1, fill(true, 4))
            ads_output2 = capture_print(AbstractWalkers.print_adsorption, s_lattice2, fill(true, 8))

            @test contains(ads_output1, "●")
            @test !contains(ads_output1, "➊")  # Shouldn't use numbered symbols for adsorption
            @test contains(ads_output2, "●")
            @test contains(ads_output2, "Layer")

            # Test multi-component print
            lattice2 = MLattice{2,SquareLattice}(
                supercell_dimensions = (2, 2, 1),
                components = [[true, false, false, false], [false, true, false, false]]
            )

            # Test print_layer_multi_comp()
            merged = AbstractWalkers.merge_components(lattice2)
            print_layer_output2 = capture_print(AbstractWalkers.print_layer, lattice2, merged)
            @test contains(print_layer_output2, "➊")
            @test contains(print_layer_output2, "➋")

            # Test occupation printing
            occ_output2 = capture_print(AbstractWalkers.print_occupation, lattice2)
            @test contains(occ_output2, "➊")
            @test contains(occ_output2, "➋")

            # Test adsorption printing
            ads_output = capture_print(AbstractWalkers.print_adsorption, lattice2, fill(true, 4))
            @test contains(ads_output, "●")
            @test !contains(ads_output, "➊")  # Shouldn't use numbered symbols for adsorption

            # Extra cases
            ## Test occupation printing for multiple layers
            lattice_3d = MLattice{1,SquareLattice}(
                supercell_dimensions = (2, 2, 2),
                components = [fill(true, 8)]
            )
            occ_output_3d = capture_print(AbstractWalkers.print_occupation, lattice_3d)
            @test contains(occ_output_3d, "Layer 1")
            @test contains(occ_output_3d, "Layer 2")
        end

        @testset "component functions utility tests" begin
            
            # Create test system with multiple species
            at = FreeBirdIO.generate_multi_type_random_starting_config(
                10.0,                              # box size
                [2,1,3,4,5,6],                    # number of particles
                particle_types=[:H,:O,:H,:Fe,:Au,:Cl]
            )

            @testset "split_components" begin
                # Test splitting by particle numbers
                components = split_components(at, [2,1,3,4,5,6])
                @test length(components) == 6
                @test length(components[1]) == 2  # H₂
                @test length(components[2]) == 1  # O
                @test length(components[3]) == 3  # H₃
                @test length(components[4]) == 4  # Fe₄
                @test length(components[5]) == 5  # Au₅
                @test length(components[6]) == 6  # Cl₆
                
                # Check bounding box and boundary conditions preserved
                @test all(cell_vectors(c) == cell_vectors(at) for c in components)
                @test all(periodicity(c) == periodicity(at) for c in components)
            end

            @testset "split_components_by_chemical_species" begin
                components = AbstractWalkers.split_components_by_chemical_species(at)
                @test length(components) == 5  # H, O, Cl, Fe, Au
                
                # Check species counts
                species_counts = map(length, components)
                @test species_counts[1] == 5   # H (all H atoms)
                @test species_counts[2] == 1   # O
                @test species_counts[3] == 6   # Cl
                @test species_counts[4] == 4   # Fe
                @test species_counts[5] == 5   # Au
                
                # Verify system properties preserved
                @test all(cell_vectors(c) == cell_vectors(at) for c in components)
                @test all(periodicity(c) == periodicity(at) for c in components)
            end

            @testset "check_num_components" begin
                # Valid case
                @test nothing === AbstractWalkers.check_num_components(3, [2,3,4], [true,false,true])
                
                # Error cases
                @test_throws ArgumentError AbstractWalkers.check_num_components(3, [2,3], [true,false,true])
                @test_throws ArgumentError AbstractWalkers.check_num_components(3, [2,3,4], [true,false])
            end

            @testset "sort_components_by_atomic_number" begin
                # Test with merge_same_species=true
                list_num_par, new_system = sort_components_by_atomic_number(at)
                @test length(list_num_par) == 5  # H, O, Cl, Fe, Au
                @test list_num_par[1] == 5       # All H atoms merged
                @test length(new_system) == length(at)
                
                # Test with merge_same_species=false
                list_num_par, new_system = sort_components_by_atomic_number(at, merge_same_species=false)
                @test length(list_num_par) == 6  # H₂, H₃, O, Cl, Fe, Au
                @test length(new_system) == length(at)
                @test cell_vectors(new_system) == cell_vectors(at)
                @test periodicity(new_system) == periodicity(at)
            end
        end
    end

    @testset "empty-configuration handling tests" begin
        box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
        pbc = (true, true, true)
        at1 = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
        empty_at = FastSystem(cell_vectors(at1), periodicity(at1),
                              empty(position(at1, :)), empty(species(at1, :)), empty(mass(at1, :)))
        at2 = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å", :Ar => [3.0, 3.0, 3.0]u"Å"], box, pbc))

        @testset "helpers on empty systems and zero-count components" begin
            comps = split_components(empty_at, [0])
            @test length(comps) == 1
            @test length(comps[1]) == 0
            @test cell_vectors(comps[1]) == cell_vectors(empty_at)
            @test periodicity(comps[1]) == periodicity(empty_at)

            mixed = split_components(at2, [0, 2])
            @test map(length, mixed) == [0, 2]
            @test position(mixed[2], :) == position(at2, :)
            @test species(mixed[2], :) == species(at2, :)

            @test AbstractWalkers.split_components_by_chemical_species(empty_at) == FastSystem[]

            list_num_par, sorted = sort_components_by_atomic_number(empty_at)
            @test isempty(list_num_par)
            @test length(sorted) == 0
            @test cell_vectors(sorted) == cell_vectors(empty_at)

            list_num_par_nm, sorted_nm = sort_components_by_atomic_number(empty_at, merge_same_species=false)
            @test isempty(list_num_par_nm)
            @test length(sorted_nm) == 0
        end

        @testset "nonempty outputs unchanged by the retyped helpers" begin
            mixed_at = FastSystem(atomic_system([:H => [1.0, 1.0, 1.0]u"Å",
                                                 :O => [2.0, 2.0, 2.0]u"Å",
                                                 :H => [3.0, 3.0, 3.0]u"Å"], box, pbc))
            comps = AbstractWalkers.split_components_by_chemical_species(mixed_at)
            @test map(length, comps) == [2, 1]  # H (Z = 1) before O (Z = 8)
            @test position(comps[1], 1) == position(mixed_at, 1)
            @test position(comps[1], 2) == position(mixed_at, 3)
            @test position(comps[2], 1) == position(mixed_at, 2)
        end

        @testset "zero-particle AtomWalker{1} round-trip" begin
            w = AtomWalker{1}(empty_at)
            @test w.list_num_par == [0]
            @test w.frozen == [false]
            lj = LJParameters()
            @test interacting_energy(w.configuration, lj) == 0.0u"eV"
            @test interacting_energy(w.configuration, lj, w.list_num_par, w.frozen) == 0.0u"eV"
            @test frozen_energy(w.configuration, lj, w.list_num_par, [true]) == 0.0u"eV"
        end

        @testset "convenience constructor rejects empty input legibly" begin
            @test_throws ArgumentError AtomWalker(empty_at)
        end

        @testset "liveset accepts a zero-particle walker" begin
            w0 = AtomWalker{1}(empty_at)
            w2 = AtomWalker{1}(at2)
            ls = LJAtomWalkers([w0, w2], LJParameters())
            @test length(ls.walkers) == 2
            @test ls.walkers[1].energy == 0.0u"eV"
            @test isfinite(ustrip(ls.walkers[2].energy))
        end
    end

    @testset "shared-geometry walkers and clones" begin
        using Random
        using Unitful
        using DataFrames

        sg_lat() = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1],
            components=[[false for _ in 1:16]], adsorptions=:full)
        sg_ham = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")

        @testset "share_geometry constructor" begin
            src = sg_lat()
            Random.seed!(98001)
            for i in 1:16
                src.components[1][i] = rand() < 0.5
            end
            occ = [copy(v) for v in src.components]
            cp = MLattice{1,SquareLattice}(Val(:share_geometry), src, occ)
            @test cp.lattice_vectors === src.lattice_vectors
            @test cp.positions === src.positions
            @test cp.basis === src.basis
            @test cp.supercell_dimensions === src.supercell_dimensions
            @test cp.periodicity === src.periodicity
            @test cp.cutoff_radii === src.cutoff_radii
            @test cp.neighbors === src.neighbors
            @test cp.adsorptions === src.adsorptions
            @test cp.components !== src.components
            @test cp.components[1] !== src.components[1]
            @test interacting_energy(cp, sg_ham) ==
                  interacting_energy(src, sg_ham)
            @test_throws ArgumentError MLattice{1,SquareLattice}(
                Val(:share_geometry), src, Vector{Bool}[])
        end

        @testset "replicate_walkers" begin
            tmpl = sg_lat()
            ws = replicate_walkers(tmpl, 5)
            @test length(ws) == 5
            @test all(w.configuration.neighbors === tmpl.neighbors
                      for w in ws)
            @test all(w.configuration.positions === tmpl.positions
                      for w in ws)
            @test length(unique(objectid(w.configuration.components[1])
                                for w in ws)) == 5
            @test all(w.energy == 0.0u"eV" && w.iter == 0 for w in ws)
            # independent occupancies: mutating one walker leaves the rest
            ws[1].configuration.components[1][1] = true
            @test !ws[2].configuration.components[1][1]
        end

        @testset "same-seed driver equivalence and aliasing safety" begin
            # A replicate_walkers live set matches the deepcopy idiom
            # digit-for-digit: both drivers re-initialize occupancies on
            # entry, so the random streams coincide
            sg_save = SaveEveryN("t_sg.csv", "t_sg.traj", "t_sg.ls",
                                 1000000, 1000000, 1000000)
            sg_cleanup() = rm.(["t_sg.csv", "t_sg.traj", "t_sg.ls"],
                               force=true)
            function sg_run(build)
                Random.seed!(98002)
                tmpl = sg_lat()
                ws = build(tmpl)
                ls = LatticeGasWalkers(ws, sg_ham; assign_energy=false)
                p = IdealGasReferencedGCNSParameters(mc_steps=20,
                    reference_fugacity=1.0, energy_perturbation=1e-9)
                d, lsx, _ = ideal_gas_referenced_nested_sampling(ls, p,
                    Int64(50), MCGrandCanonicalMoves(), sg_save)
                sg_cleanup()
                return d, lsx
            end
            dS, lsS = sg_run(t -> replicate_walkers(t, 10))
            dD, lsD = sg_run(t -> [LatticeWalker(deepcopy(t),
                                                 energy=0.0u"eV", iter=0)
                                   for _ in 1:10])
            @test dS.iter == dD.iter
            @test dS.emax == dD.emax
            @test dS.num_particles == dD.num_particles
            @test [w.energy.val for w in lsS.walkers] ==
                  [w.energy.val for w in lsD.walkers]
            # aliasing safety: after the run every live walker still shares
            # the one geometry (nothing rebinds or copies it back)
            nb1 = lsS.walkers[1].configuration.neighbors
            @test all(w.configuration.neighbors === nb1
                      for w in lsS.walkers)
        end

        @testset "clone allocation guard" begin
            # Compile-warmed shared clone at M = 4096 under a generous fixed
            # byte ceiling (never time-based): the payload is one occupancy
            # vector plus constant-size shells
            big = MLattice{1,SquareLattice}(lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(64, 64, 1),
                periodicity=(true, true, false), cutoff_radii=[1.1],
                components=[[false for _ in 1:4096]], adsorptions=:full)
            wk = LatticeWalker(big, energy=0.0u"eV", iter=0)
            FreeBird.SamplingSchemes._clone_walker_shared_geometry(wk)
            bytes = @allocated FreeBird.SamplingSchemes._clone_walker_shared_geometry(wk)
            @test bytes < 50_000
        end
    end

end
