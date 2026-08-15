# Exactness coverage for the O(z) single-site flip deltas and the
# supports_site_deltas trait. All assertions are exact or fixed-tolerance;
# there are no statistical gates. The atol 1e-13 eV carries a 12x margin
# over the measured worst case of the prototype (<= 8e-15 eV over 200
# random flips per size).

# Minimal trait fixture: testset bodies cannot define structs, so it sits
# at the file's top level (identical re-definition on a double include is
# safe).
struct SFDMinimalHam <: FreeBird.AbstractHamiltonians.ClassicalHamiltonian end

@testset "Lattice site-flip deltas" begin
    using Random
    using Unitful

    function sfd_maxdev(lat, h, seed; nflips=200)
        Random.seed!(seed)
        maxdev = 0.0
        M = length(lat.components[1])
        for _ in 1:nflips
            s = rand(1:M)
            e0 = interacting_energy(lat, h)
            d = site_flip_delta(lat, h, s)
            lat.components[1][s] = !lat.components[1][s]
            e1 = interacting_energy(lat, h)
            lat.components[1][s] = !lat.components[1][s]
            maxdev = max(maxdev, abs(ustrip(u"eV", (e1 - e0) - d)))
        end
        return maxdev
    end

    function sfd_square(L; cutoffs=[1.1], adsorptions=:full,
                        image_multiplicity=false)
        MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false), cutoff_radii=cutoffs,
            components=[[false for _ in 1:L*L]], adsorptions=adsorptions,
            image_multiplicity=image_multiplicity)
    end

    function sfd_fill!(lat, seed, theta)
        Random.seed!(seed)
        for i in eachindex(lat.components[1])
            lat.components[1][i] = rand() < theta
        end
        return lat
    end

    h1 = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
    h2 = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

    @testset "exactness against full-energy differences" begin
        # square, one shell, full adsorption
        @test sfd_maxdev(sfd_fill!(sfd_square(4), 92_001, 0.5), h1,
                         93_001) <= 1e-13
        # square, two shells
        @test sfd_maxdev(sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]),
                                   92_002, 0.5), h2, 93_002) <= 1e-13
        # partial adsorption mask
        @test sfd_maxdev(sfd_fill!(sfd_square(4;
                                              adsorptions=[1, 2, 3, 5, 8, 13]),
                                   92_003, 0.3), h1, 93_003) <= 1e-13
        # triangular geometry
        tri = SLattice{TriangularLattice}(supercell_dimensions=(6, 4, 1),
            cutoff_radii=[1.1], components=[[1, 5, 10, 20, 30, 40]],
            adsorptions=:full)
        @test sfd_maxdev(sfd_fill!(tri, 92_004, 0.5), h1, 93_004) <= 1e-13
        # image-multiplicity small cell: self-image entries and duplicated
        # image bonds sit under the delta's j == site convention
        imcell = sfd_square(2; image_multiplicity=true)
        @test sfd_maxdev(sfd_fill!(imcell, 92_005, 0.5), h1, 93_005) <= 1e-13
        # site-field wrapper
        fld = collect(0.001 .* (1:16)) .* u"eV"
        hsf = SiteFieldLatticeHamiltonian(h1, fld)
        @test sfd_maxdev(sfd_fill!(sfd_square(4), 92_006, 0.5), hsf,
                         93_006) <= 1e-13
        # multi-component Hamiltonian on a single-component lattice
        mlham = MLatticeHamiltonian(1, [h1])
        @test sfd_maxdev(sfd_fill!(sfd_square(4), 92_007, 0.5), mlham,
                         93_007) <= 1e-13
    end

    @testset "sign symmetry is exact" begin
        # The accumulator visits the same terms with the opposite sign, so
        # the back-flip delta is the exact floating-point negation
        lat = sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]), 92_008, 0.5)
        Random.seed!(93_008)
        for _ in 1:50
            s = rand(1:36)
            d_fwd = site_flip_delta(lat, h2, s)
            lat.components[1][s] = !lat.components[1][s]
            d_back = site_flip_delta(lat, h2, s)
            lat.components[1][s] = !lat.components[1][s]
            @test d_back == -d_fwd
        end
    end

    @testset "swap composition" begin
        # Two sequentially composed flips, the second evaluated on the
        # intermediate configuration, reproduce the full-energy difference
        # of the swap at the single-flip rounding class, including adjacent
        # origin-destination pairs
        lat = sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]), 92_009, 0.5)
        Random.seed!(93_009)
        maxdev = 0.0
        for trial in 1:100
            i = rand(1:36)
            # force adjacency on odd trials: destination = a first-shell
            # neighbor entry of the origin
            j = isodd(trial) ? lat.neighbors[i][1][1] : rand(1:36)
            e0 = interacting_energy(lat, h2)
            d = site_flip_delta(lat, h2, i)
            lat.components[1][i] = !lat.components[1][i]
            d += site_flip_delta(lat, h2, j)
            lat.components[1][j] = !lat.components[1][j]
            e1 = interacting_energy(lat, h2)
            # restore
            lat.components[1][j] = !lat.components[1][j]
            lat.components[1][i] = !lat.components[1][i]
            maxdev = max(maxdev, abs(ustrip(u"eV", (e1 - e0) - d)))
        end
        @test maxdev <= 1e-13
    end

    @testset "trait contract and shell guard" begin
        fld = collect(0.001 .* (1:16)) .* u"eV"
        clham = ClusterLatticeHamiltonian(h1,
            [ClusterInteraction(0.1u"eV", [(1, 2, 3)])])
        @test supports_site_deltas(h1)
        @test supports_site_deltas(MLatticeHamiltonian(1, [h1]))
        @test supports_site_deltas(SiteFieldLatticeHamiltonian(h1, fld))
        @test !supports_site_deltas(SFDMinimalHam())
        @test !supports_site_deltas(clham)
        # the wrapper delegates to its base
        @test !supports_site_deltas(SiteFieldLatticeHamiltonian(clham, fld))
        # the delta enforces the same shell-count rule as the full sweep
        lat1 = sfd_square(4)
        @test_throws ArgumentError site_flip_delta(lat1, h2, 1)
    end
end
