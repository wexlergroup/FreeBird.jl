@testset "exact_enumeration.jl tests" begin
    @testset "unique_permutations function tests" begin
        # Test empty input
        @test SamplingSchemes.unique_permutations(Int[]) == Vector{Int}[]
        
        # Test single element
        @test SamplingSchemes.unique_permutations([1]) == [[1]]
        
        # Test two unique elements
        @test sort(SamplingSchemes.unique_permutations([1,2])) == sort([[1,2], [2,1]])
        
        # Test two identical elements
        @test SamplingSchemes.unique_permutations([1,1]) == [[1,1]]
        
        # Test three elements with duplicates
        result3 = SamplingSchemes.unique_permutations([1,1,2])
        @test length(result3) == 3
        @test sort(result3) == sort([[1,1,2], [1,2,1], [2,1,1]])
        
        # Test string input
        @test sort(SamplingSchemes.unique_permutations(['a','b'])) == sort([['a','b'], ['b','a']])
        
        # Test longer sequence with multiple duplicates
        result4 = SamplingSchemes.unique_permutations([1,1,2,2])
        @test length(result4) == 6
        @test sort(result4) == sort([
            [1,1,2,2], [1,2,1,2], [1,2,2,1],
            [2,1,1,2], [2,1,2,1], [2,2,1,1]
        ])
        
        # Test type stability
        @test eltype(SamplingSchemes.unique_permutations([1,2])) == Vector{Int}
        @test eltype(SamplingSchemes.unique_permutations(['a','b'])) == Vector{Char}
    end

    @testset "enumerate_lattices function tests" begin
        @testset "MLattice as init_lattice cases" begin
            # 2x1 square lattice with 2 components
            lattice = MLattice{2,SquareLattice}(
                supercell_dimensions=(2, 1, 1),
                components=[[1], [2]]
            )
            results = SamplingSchemes.enumerate_lattices(lattice)
            
            @test length(results) == 2
            @test all(r isa MLattice{2,SquareLattice} for r in results)
            @test all(r.lattice_vectors == lattice.lattice_vectors for r in results)
            @test all(r.supercell_dimensions == (2,1,1) for r in results)
            @test all(sum(r.components[1]) + sum(r.components[2]) == 2 for r in results)
            
            # 2x2 triangular lattice with 3 components
            lattice = MLattice{3,TriangularLattice}(
                supercell_dimensions=(2, 2, 1),
                components=[[1], [2], [3]]
            )
            results = SamplingSchemes.enumerate_lattices(lattice)
            
            @test length(results) == 336  # (8! / (3! * 5!)) * 3 * 2 * 1
            @test all(r isa MLattice{3,TriangularLattice} for r in results)
            @test all(r.basis == lattice.basis for r in results)
            @test all(sum(r.components[1]) == 1 for r in results)
            @test all(sum(r.components[2]) == 1 for r in results)
            @test all(sum(r.components[3]) == 1 for r in results)
        end

        @testset "SLattice as init_lattice cases" begin
            # 2x1 square lattice with 1 occupied site
            lattice = SLattice{SquareLattice}(
                supercell_dimensions=(2, 1, 1),
                components=[[true, false]]
            )
            results = SamplingSchemes.enumerate_lattices(lattice)
            
            @test length(results) == 2  # C(2,1)
            @test all(r isa SLattice{SquareLattice} for r in results)
            @test all(sum(r.components[1]) == 1 for r in results)
            
            # 2x2 triangular lattice with 2 occupied sites
            lattice = SLattice{TriangularLattice}(
                supercell_dimensions=(2, 2, 1),
                components=[[true, true, true, true, 
                false, false, false, false]]
            )
            results = SamplingSchemes.enumerate_lattices(lattice)
            
            @test length(results) == 70  # C(8,4)
            @test all(r isa SLattice{TriangularLattice} for r in results)
            @test all(sum(r.components[1]) == 4 for r in results)

            # Verify thread safety
            results_set = Set(map(r -> findall(r.components[1]), results))
            @test length(results_set) == 70
        end
    end


    @testset "exact_enumeration function tests" begin
        # Test 2-component system
        on_site = [-0.04, -0.03]u"eV"
        neighbor_ints = [[-0.01, -0.0025], [-0.015, -0.002], [-0.02, -0.003]].*u"eV"
        hams = [GenericLatticeHamiltonian(on_site[i], neighbor_ints[i]) for i in 1:2]
        push!(hams, GenericLatticeHamiltonian(on_site[2], neighbor_ints[3]))
        h = MLatticeHamiltonian(2, hams)
        
        lattice = MLattice{2,SquareLattice}(
            supercell_dimensions=(2, 2, 1),
            components=[[1], [4]]
        )
        
        df, ls = exact_enumeration(lattice, h)
        
        @test size(df, 1) == 12
        @test all(x -> unit(x) == u"eV", df.energy)
        @test all(c -> length(c) == 2, df.config)
        
        # Test 3-component system
        on_site = [-0.04, -0.03, -0.02]u"eV"
        neighbor_ints = [
            [-0.01, -0.0025], 
            [-0.015, -0.002], 
            [-0.02, -0.003],
            [-0.012, -0.0022],
            [-0.018, -0.0028],
            [-0.016, -0.0024]
        ].*u"eV"
        
        hams = [GenericLatticeHamiltonian(on_site[i], neighbor_ints[i]) for i in 1:3]
        for i in 4:6
            push!(hams, GenericLatticeHamiltonian(on_site[2], neighbor_ints[i]))
        end
        h = MLatticeHamiltonian(3, hams)
        
        lattice = MLattice{3,SquareLattice}(
            supercell_dimensions=(2, 2, 1),
            components=[[true, false, false, false], 
                    [false, true, true, false],
                    [false, false, false, true]]
        )
        
        df, ls = exact_enumeration(lattice, h)
        
        @test size(df, 1) == 12
        @test all(x -> unit(x) == u"eV", df.energy)
        @test all(sum(c[1]) == 1 for c in df.config)
        @test all(sum(c[2]) == 2 for c in df.config)
        @test all(sum(c[3]) == 1 for c in df.config)
        
        # Test energy differences
        @test length(unique(df.energy)) > 1
    end

    @testset "grand_canonical_exact_enumeration tests" begin
        L = 3
        M = L * L
        lat = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:M]],
            adsorptions=:full
        )
        ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

        E, N = grand_canonical_exact_enumeration(lat, ham)

        # Every occupation, once, in mask order.
        @test length(E) == 2^M
        @test length(N) == 2^M
        @test all(N[i] == count_ones(i - 1) for i in eachindex(N))
        @test minimum(N) == 0
        @test maximum(N) == M
        @test count(==(0), N) == 1          # the empty lattice
        @test count(==(M), N) == 1          # the full one

        # Bit-identical to the hand-rolled 2^M loop this replaces in the
        # grand-canonical test files. Exact equality, not ≈: it is the same
        # deterministic computation, and anything else would mean the strided
        # threading had perturbed it.
        ref_E = similar(E)
        ref_N = zeros(Int, 2^M)
        tmp = deepcopy(lat)
        for mask in 0:(2^M - 1)
            for site in 1:M
                tmp.components[1][site] = ((mask >> (site - 1)) & 1) == 1
            end
            ref_E[mask + 1] = interacting_energy(tmp, ham)
            ref_N[mask + 1] = sum(tmp.components[1])
        end
        @test E == ref_E
        @test N == ref_N

        # The caller's lattice is a template, not a scratch buffer.
        @test all(.!lat.components[1])

        # A grand-canonical average assembled from the spectrum: μ and β belong
        # to the ensemble, not to the enumeration.
        μ, β = -0.05, 1.0 / (8.617333262e-5 * 300.0)
        Ω = ustrip.(E) .- μ .* N
        w = exp.(-β .* (Ω .- minimum(Ω)))
        @test 0.0 <= sum(w .* N) / sum(w) <= M

        # Guard against asking for 2^40 states by accident.
        @test_throws ArgumentError grand_canonical_exact_enumeration(lat, ham; max_sites=4)
    end

end