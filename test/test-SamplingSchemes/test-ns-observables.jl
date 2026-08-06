@testset "NS observables" begin
    using Random

    # Shared single-component square-lattice builder
    function obs_square_lattice(d1, d2)
        MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(d1, d2, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:d1*d2]],
            adsorptions=:full)
    end

    # Shared single-component triangular-lattice builder (NN shell only)
    function obs_tri_lattice(d1, d2)
        MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            supercell_dimensions=(d1, d2, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:2*d1*d2]],
            adsorptions=:full)
    end

    # Shared single-component square slab builder, d3 layers of d1 x d2
    # sites. At d1 = d2 = 2 the first shell wraps the periodic cell (each
    # in-plane pair is connected through both images), which the default
    # minimum-image lists warn about at construction; build with
    # image_multiplicity=true so construction stays silent. The doubled
    # in-plane bonds are energetically irrelevant everywhere the builder
    # is used (J = 0 throughout, and the layer helpers never read
    # neighbor lists).
    function obs_slab_lattice(d1, d2, d3)
        MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(d1, d2, d3),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:d1*d2*d3]],
            adsorptions=:full,
            image_multiplicity=true)
    end

    # ================================================================
    @testset "order_parameter_c2x2" begin
        lat = obs_square_lattice(4, 4)

        # Empty and full lattices carry no c(2x2) order
        @test order_parameter_c2x2(lat) == 0.0
        lat.components[1] .= true
        @test order_parameter_c2x2(lat) == 0.0

        # Single particle: |±1| / M
        lat.components[1] .= false
        lat.components[1][1] = true
        @test order_parameter_c2x2(lat) == 1 / 16

        # Perfect checkerboards (both sublattices) at half filling: 1/2.
        # Site ordering: dimension 1 fastest, so i0 = (s-1) % d1, j0 = (s-1) ÷ d1.
        checker = [iseven(((s - 1) % 4) + ((s - 1) ÷ 4)) for s in 1:16]
        lat.components[1] .= checker
        @test order_parameter_c2x2(lat) == 0.5
        lat.components[1] .= .!checker
        @test order_parameter_c2x2(lat) == 0.5

        # A column stripe alternates parity along the column and cancels
        lat.components[1] .= [(s - 1) % 4 == 0 for s in 1:16]
        @test order_parameter_c2x2(lat) == 0.0

        # Odd in-plane dimensions: the sublattices do not tile evenly
        @test_throws ArgumentError order_parameter_c2x2(obs_square_lattice(3, 4))
        @test_throws ArgumentError order_parameter_c2x2(obs_square_lattice(4, 3))

        # Three-dimensional supercell
        lat3d = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 2),
            periodicity=(true, true, true),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:32]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_c2x2(lat3d)

        # Multi-site basis
        lat2b = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)],
            supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[0.8],
            components=[[false for _ in 1:32]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_c2x2(lat2b)
    end

    # ================================================================
    obs_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
    obs_save = SaveEveryN("t_obs.csv", "t_obs.traj", "t_obs.ls",
                          1000000, 1000000, 1000000)
    obs_cleanup() = rm.(["t_obs.csv", "t_obs.traj", "t_obs.ls"], force=true)
    count_N(cfg) = Float64(sum(cfg.components[1]))

    # ================================================================
    @testset "order_parameter_sqrt3" begin
        # Independent geometric decoder: recover the triangular integer
        # coordinates (m, n) by inverting [t1 t2] on the stored positions and
        # color with c = (m - n) mod 3, so the labels are derived from the
        # geometry rather than from the index arithmetic under test
        function sqrt3_labels(lat)
            a = lat.lattice_vectors[1, 1]
            map(1:size(lat.positions, 1)) do s
                x, y = lat.positions[s, 1], lat.positions[s, 2]
                n = round(Int, 2 * y / (sqrt(3) * a))
                m = round(Int, x / a - n / 2)
                mod(m - n, 3)
            end
        end

        for (d1, d2) in [(3, 3), (6, 2)]
            lat = obs_tri_lattice(d1, d2)
            M = 2 * d1 * d2
            labels = sqrt3_labels(lat)
            # The hardcoded decoder inside order_parameter_sqrt3, recomputed
            # here, must match the geometric labels site by site
            hard = [begin
                        b = (s - 1) % 2
                        ci = ((s - 1) ÷ 2) % d1
                        b == 0 ? ci % 3 : (ci + 2) % 3
                    end for s in 1:M]
            @test labels == hard
            # Each sublattice holds exactly M/3 sites
            @test [count(==(c), labels) for c in 0:2] == fill(M ÷ 3, 3)
            # A proper 3-coloring: no first-shell neighbor pair shares a
            # label (this also pins the neighbor lists)
            @test all(labels[nb] != labels[s]
                      for s in 1:M for nb in lat.neighbors[s][1])
        end

        # Exact values on the 18-site commensurate cell
        lat = obs_tri_lattice(3, 3)
        labels = sqrt3_labels(lat)

        # Empty and full lattices carry no sqrt3 order (1 + ω + ω² = 0)
        @test order_parameter_sqrt3(lat) == 0.0
        lat.components[1] .= true
        @test order_parameter_sqrt3(lat) == 0.0

        # Each pure sublattice state, built from the geometric labels: 1/3
        for c in 0:2
            lat.components[1] .= (labels .== c)
            @test order_parameter_sqrt3(lat) == 1 / 3
        end

        # Single particle: 1/M
        lat.components[1] .= false
        lat.components[1][1] = true
        @test order_parameter_sqrt3(lat) == 1 / 18

        # Mixed state with sublattice counts (2, 1, 0): |2 + ω| / 18 = √3/18
        lat.components[1] .= false
        a_sites = findall(==(0), labels)
        b_sites = findall(==(1), labels)
        lat.components[1][a_sites[1]] = true
        lat.components[1][a_sites[2]] = true
        lat.components[1][b_sites[1]] = true
        @test order_parameter_sqrt3(lat) ≈ sqrt(3) / 18 rtol = 1e-12

        # The shipped default (4, 2, 1) cell is incommensurate: d1 % 3 != 0
        @test_throws ArgumentError order_parameter_sqrt3(obs_tri_lattice(4, 2))

        # Three-dimensional supercell
        lat3d = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            supercell_dimensions=(3, 3, 2),
            periodicity=(true, true, true),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:36]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_sqrt3(lat3d)

        # One-site basis
        lat1b = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(3, 3, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:9]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_sqrt3(lat1b)

        # A distorted basis breaks the tripartition arithmetic
        latdb = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.3, 0.3, 0.0)],
            supercell_dimensions=(3, 3, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:18]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_sqrt3(latdb)

        # lattice_constant != 1 scales lattice_vectors but not the default
        # basis, so the geometry is no longer triangular; the basis-values
        # guard rejects it
        latlc = MLattice{1,TriangularLattice}(
            lattice_constant=2.0,
            supercell_dimensions=(3, 3, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.2],
            components=[[false for _ in 1:18]],
            adsorptions=:full)
        @test_throws ArgumentError order_parameter_sqrt3(latlc)

        # A short seeded canonical NS run records psi and its square through
        # the moment-callback pattern from the docstring
        Random.seed!(21)
        tri_ham = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
        walkers = [LatticeWalker(deepcopy(obs_tri_lattice(3, 3)),
                                 energy=0.0u"eV", iter=0) for _ in 1:20]
        # Fixed-N ladder at N = 6: place 6 particles per walker
        for w in walkers
            occ = vcat(fill(true, 6), fill(false, 12))
            shuffle!(occ)
            w.configuration.components[1] .= occ
        end
        ls = LatticeGasWalkers(walkers, tri_ham; perturb_energy=1e-9)
        params = LatticeNestedSamplingParameters(mc_steps=30,
            energy_perturbation=1e-9, allowed_fail_count=100000)

        df, _, _ = nested_sampling(ls, params, Int64(200),
            MCRandomWalkClone(), obs_save;
            observables=[:psi => order_parameter_sqrt3,
                         :psi2 => (cfg -> order_parameter_sqrt3(cfg)^2)])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "psi", "psi2"]
        @test all(0.0 .<= df.psi .<= 1 / 3)
        @test all(isapprox.(df.psi2, df.psi .^ 2; rtol=1e-12))

        # Schema is unchanged when observables are not requested
        Random.seed!(21)
        walkers2 = [LatticeWalker(deepcopy(obs_tri_lattice(3, 3)),
                                  energy=0.0u"eV", iter=0) for _ in 1:20]
        for w in walkers2
            occ = vcat(fill(true, 6), fill(false, 12))
            shuffle!(occ)
            w.configuration.components[1] .= occ
        end
        ls2 = LatticeGasWalkers(walkers2, tri_ham; perturb_energy=1e-9)
        params2 = LatticeNestedSamplingParameters(mc_steps=30,
            energy_perturbation=1e-9, allowed_fail_count=100000)
        df2, _, _ = nested_sampling(ls2, params2, Int64(50),
            MCRandomWalkClone(), obs_save)
        obs_cleanup()
        @test names(df2) == ["iter", "emax"]
    end

    # ================================================================
    @testset "bragg_amplitude and order_parameter_stripe" begin
        lat = obs_square_lattice(4, 4)
        # Occupation from a predicate on the integer coordinates; site
        # ordering: dimension 1 fastest, so x = (s-1) % d1, y = (s-1) ÷ d1
        function fill_xy!(l, pred)
            d1, d2 = l.supercell_dimensions[1], l.supercell_dimensions[2]
            l.components[1] .= [pred((s - 1) % d1, ((s - 1) ÷ d1) % d2)
                                for s in 1:d1*d2]
            l
        end

        # k = (π, π) is order_parameter_c2x2: the two code paths agree on
        # random occupations, and indices wrap modulo the dimensions
        Random.seed!(1234)
        for _ in 1:20
            lat.components[1] .= rand(Bool, 16)
            @test bragg_amplitude(lat, 2, 2) == order_parameter_c2x2(lat)
            @test bragg_amplitude(lat, 1, 3) == bragg_amplitude(lat, 5, 3)
            @test bragg_amplitude(lat, 1, 3) == bragg_amplitude(lat, -3, 3)
        end

        # Perfect (2x1) stripes, both orientations, both phases: 1/2, with
        # order_parameter_c2x2 blind to all four (the :49 anti-case, from
        # the other side); the perfect checkerboard is the converse
        for pred in [(x, y) -> x % 2 == 0, (x, y) -> x % 2 == 1,
                     (x, y) -> y % 2 == 0, (x, y) -> y % 2 == 1]
            fill_xy!(lat, pred)
            @test order_parameter_stripe(lat) ≈ 1 / 2 atol = 1e-12
            @test order_parameter_c2x2(lat) == 0.0
        end
        fill_xy!(lat, (x, y) -> iseven(x + y))
        @test order_parameter_stripe(lat) == 0.0

        # Period-4 axial double stripes, both orientations, all four phases,
        # on 4x4 and 8x8: sqrt(2)/4 at the stripe harmonic, exactly 0 at the
        # period-2 harmonic (wrong-harmonic blindness)
        for (d1, d2) in [(4, 4), (8, 8)]
            latp = obs_square_lattice(d1, d2)
            for p in 0:3, pred in [(x, y) -> mod(x - p, 4) < 2,
                                   (x, y) -> mod(y - p, 4) < 2]
                fill_xy!(latp, pred)
                @test order_parameter_stripe(latp, period=4) ≈ sqrt(2) / 4 atol = 1e-12
                @test order_parameter_stripe(latp, period=2) == 0.0
            end
        end

        # Block and diagonal period-4 states: the axial stripe order
        # parameter is exactly 0 while the documented diagonal quadrature
        # resolves the k = (π/2, ±π/2) manifold at sqrt(2)/4
        diag_quad(cfg) = sqrt(bragg_amplitude(cfg, 1, 1)^2 +
                              bragg_amplitude(cfg, 1, -1)^2)
        fill_xy!(lat, (x, y) -> x ÷ 2 == y ÷ 2)
        @test order_parameter_stripe(lat, period=4) == 0.0
        @test bragg_amplitude(lat, 1, 1) ≈ 1 / 4 atol = 1e-12
        @test bragg_amplitude(lat, 1, -1) ≈ 1 / 4 atol = 1e-12
        @test diag_quad(lat) ≈ sqrt(2) / 4 atol = 1e-12
        fill_xy!(lat, (x, y) -> mod(x + y, 4) < 2)
        @test order_parameter_stripe(lat, period=4) == 0.0
        @test bragg_amplitude(lat, 1, 1) ≈ sqrt(2) / 4 atol = 1e-12
        @test bragg_amplitude(lat, 1, -1) == 0.0
        @test diag_quad(lat) ≈ sqrt(2) / 4 atol = 1e-12

        # Empty and full lattices scatter nothing off the (0, 0) rod
        lat.components[1] .= false
        @test bragg_amplitude(lat, 1, 0) == 0.0
        @test order_parameter_stripe(lat) == 0.0
        lat.components[1] .= true
        for (m, n) in [(1, 0), (0, 1), (2, 2), (3, 1)]
            @test bragg_amplitude(lat, m, n) == 0.0
        end
        @test order_parameter_stripe(lat) == 0.0
        @test order_parameter_stripe(lat, period=4) == 0.0

        # Single particle: 1/M at any k, sqrt(2)/M in quadrature
        lat.components[1] .= false
        lat.components[1][1] = true
        @test bragg_amplitude(lat, 1, 0) == 1 / 16
        @test bragg_amplitude(lat, 2, 3) == 1 / 16
        @test order_parameter_stripe(lat) ≈ sqrt(2) / 16 atol = 1e-12

        # Two-site basis
        lat2b = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)],
            supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[0.8],
            components=[[false for _ in 1:32]],
            adsorptions=:full)
        @test_throws ArgumentError bragg_amplitude(lat2b, 1, 0)
        @test_throws ArgumentError order_parameter_stripe(lat2b)

        # Three-dimensional supercell
        lat3d = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 2),
            periodicity=(true, true, true),
            cutoff_radii=[1.1, 1.5],
            components=[[false for _ in 1:32]],
            adsorptions=:full)
        @test_throws ArgumentError bragg_amplitude(lat3d, 1, 0)
        @test_throws ArgumentError order_parameter_stripe(lat3d)

        # Incommensurate and degenerate periods
        @test_throws ArgumentError order_parameter_stripe(
            obs_square_lattice(4, 4), period=3)
        @test_throws ArgumentError order_parameter_stripe(
            obs_square_lattice(6, 4), period=4)
        @test_throws ArgumentError order_parameter_stripe(
            obs_square_lattice(4, 4), period=1)
        @test_throws ArgumentError order_parameter_stripe(
            obs_square_lattice(4, 4), period=0)

        # A short seeded canonical NS run records the stripe Binder moments
        # through the moment-callback pattern from the docstring
        Random.seed!(23)
        walkers = [LatticeWalker(deepcopy(obs_square_lattice(4, 4)),
                                 energy=0.0u"eV", iter=0) for _ in 1:20]
        # Fixed-N ladder at N = 6: place 6 particles per walker
        for w in walkers
            occ = vcat(fill(true, 6), fill(false, 10))
            shuffle!(occ)
            w.configuration.components[1] .= occ
        end
        ls = LatticeGasWalkers(walkers, obs_ham; perturb_energy=1e-9)
        params = LatticeNestedSamplingParameters(mc_steps=30,
            energy_perturbation=1e-9, allowed_fail_count=100000)

        df, _, _ = nested_sampling(ls, params, Int64(200),
            MCRandomWalkClone(), obs_save;
            observables=[:stripe => order_parameter_stripe,
                         :stripe2 => (cfg -> order_parameter_stripe(cfg)^2),
                         :stripe4 => (cfg -> order_parameter_stripe(cfg)^4)])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "stripe", "stripe2", "stripe4"]
        @test all(0.0 .<= df.stripe .<= 1 / 2)
        @test all(isapprox.(df.stripe2, df.stripe .^ 2; rtol=1e-12))
        @test all(isapprox.(df.stripe4, df.stripe .^ 4; rtol=1e-12))
    end

    # ================================================================
    @testset "triangular bragg_amplitude and order_parameter_p2x2" begin
        # Axial coordinates from the stored positions: r = i·t₁ + j·t₂ with
        # t₁ = (1, 0)·a and t₂ = (1/2, √3/2)·a, so j = 2y/(√3·a) and
        # i = x/a − j/2. Occupation predicates on (i, j) build the ordered
        # phases geometrically, independent of the index arithmetic under
        # test: p(2×2) = {i, j both even}, a p(2×1) row phase = {j even}
        # (plus its two rotations), the √3×√3 state = {(i − j) mod 3 == 0}.
        function axial_ij(lat)
            a = lat.lattice_vectors[1, 1]
            map(1:size(lat.positions, 1)) do s
                x, y = lat.positions[s, 1], lat.positions[s, 2]
                j = round(Int, 2 * y / (sqrt(3) * a))
                i = round(Int, x / a - j / 2)
                (i, j)
            end
        end
        fill_ij!(lat, pred) = (lat.components[1] .=
            [pred(i, j) for (i, j) in axial_ij(lat)]; lat)

        # Independent oracle: the raw phasor sum over the stored positions,
        # sharing no code with the integer-table implementation
        function tri_oracle(lat, m, n)
            a = lat.lattice_vectors[1, 1]
            d1, d2 = lat.supercell_dimensions[1], lat.supercell_dimensions[2]
            kx = 2 * pi * m / (d1 * a)
            ky = 2 * pi * n / (sqrt(3) * d2 * a)
            z = sum(exp(im * (kx * lat.positions[s, 1] + ky * lat.positions[s, 2]))
                    for s in 1:size(lat.positions, 1) if lat.components[1][s];
                    init=0.0 + 0.0im)
            return abs(z) / (2 * d1 * d2)
        end

        # The three M points on an even-dimension cell
        mpoints(d1, d2) = [(d1 ÷ 2, -(d2 ÷ 2)), (0, d2), (d1 ÷ 2, d2 ÷ 2)]

        # Oracle agreement on random occupations at the M points and at
        # generic indices, and the primitive-reciprocal wrapping identities,
        # exact by the joint mod reduction
        lat = obs_tri_lattice(6, 4)
        Random.seed!(4164)
        for _ in 1:20
            lat.components[1] .= rand(Bool, 48)
            for (m, n) in vcat(mpoints(6, 4), [(1, 0), (0, 1), (2, 3), (-1, 2)])
                @test isapprox(bragg_amplitude(lat, m, n), tri_oracle(lat, m, n);
                               atol=1e-12)
                @test bragg_amplitude(lat, m, n) == bragg_amplitude(lat, m + 6, n - 4)
                @test bragg_amplitude(lat, m, n) == bragg_amplitude(lat, m, n + 8)
            end
        end

        # Perfect-phase table on (6, 4) and (12, 6): p(2×2) puts exactly 1/4
        # on every M point (quadrature √3/4); each p(2×1) orientation puts
        # exactly 1/2 on one M point and 0 on the other two (quadrature
        # 1/2); the √3×√3 state scatters off the K points only, and the two
        # order parameters are blind to each other's phase
        for (d1, d2) in [(6, 4), (12, 6)]
            latp = obs_tri_lattice(d1, d2)
            M = 2 * d1 * d2

            fill_ij!(latp, (i, j) -> iseven(i) && iseven(j))
            @test count(latp.components[1]) == M ÷ 4
            for (m, n) in mpoints(d1, d2)
                @test bragg_amplitude(latp, m, n) == 1 / 4
            end
            @test order_parameter_p2x2(latp) ≈ sqrt(3) / 4 atol = 1e-12
            @test order_parameter_sqrt3(latp) == 0.0

            for pred in [(i, j) -> iseven(j), (i, j) -> iseven(i),
                         (i, j) -> iseven(i - j)]
                fill_ij!(latp, pred)
                @test count(latp.components[1]) == M ÷ 2
                amps = sort([bragg_amplitude(latp, m, n)
                             for (m, n) in mpoints(d1, d2)])
                @test amps[1] == 0.0 && amps[2] == 0.0 && amps[3] == 1 / 2
                @test order_parameter_p2x2(latp) ≈ 1 / 2 atol = 1e-12
                @test order_parameter_sqrt3(latp) == 0.0
            end

            fill_ij!(latp, (i, j) -> mod(i - j, 3) == 0)
            @test count(latp.components[1]) == M ÷ 3
            for (m, n) in mpoints(d1, d2)
                @test bragg_amplitude(latp, m, n) == 0.0
            end
            @test order_parameter_p2x2(latp) == 0.0
            @test order_parameter_sqrt3(latp) == 1 / 3
        end

        # Empty and full lattices scatter nothing off the reciprocal-lattice
        # rods: exactly at the M points (half-turn phases are exact), to
        # roundoff at generic indices (the phase rationals are not
        # representable); a single particle gives 1/M
        lat.components[1] .= false
        for (m, n) in vcat(mpoints(6, 4), [(1, 1)])
            @test bragg_amplitude(lat, m, n) == 0.0
        end
        @test order_parameter_p2x2(lat) == 0.0
        lat.components[1] .= true
        for (m, n) in mpoints(6, 4)
            @test bragg_amplitude(lat, m, n) == 0.0
        end
        @test bragg_amplitude(lat, 1, 1) < 1e-15
        @test order_parameter_p2x2(lat) == 0.0
        lat.components[1] .= false
        lat.components[1][1] = true
        @test bragg_amplitude(lat, 1, 1) == 1 / 48
        @test bragg_amplitude(lat, 0, 4) == 1 / 48
        @test order_parameter_p2x2(lat) ≈ sqrt(3) / 48 atol = 1e-12

        # Dispatch separation: the square methods are untouched by the new
        # triangular method
        sq = obs_square_lattice(4, 4)
        Random.seed!(4165)
        sq.components[1] .= rand(Bool, 16)
        @test bragg_amplitude(sq, 2, 2) == order_parameter_c2x2(sq)
        @test order_parameter_stripe(sq) == sqrt(bragg_amplitude(sq, 2, 0)^2 +
                                                 bragg_amplitude(sq, 0, 2)^2)

        # Guard paths: a three-dimensional supercell rejects both functions;
        # odd in-plane dimensions reject the quadrature while the
        # always-representable M₂ amplitude still evaluates; a nonstandard
        # two-site basis rejects the amplitude
        lat3d = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            supercell_dimensions=(6, 4, 2),
            periodicity=(true, true, true),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:96]],
            adsorptions=:full)
        @test_throws ArgumentError bragg_amplitude(lat3d, 1, 0)
        @test_throws ArgumentError order_parameter_p2x2(lat3d)
        for (d1, d2) in [(3, 3), (5, 4)]
            lodd = obs_tri_lattice(d1, d2)
            @test_throws ArgumentError order_parameter_p2x2(lodd)
            lodd.components[1][1] = true
            @test bragg_amplitude(lodd, 0, d2) == 1 / (2 * d1 * d2)
        end
        latnb = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)],
            supercell_dimensions=(6, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[0.75],
            components=[[false for _ in 1:48]],
            adsorptions=:full,
            image_multiplicity=true)
        @test_throws ArgumentError bragg_amplitude(latnb, 1, 0)
        @test_throws ArgumentError order_parameter_p2x2(latnb)

        # The shipped default (4, 2, 1) cell passes this function's guards
        # while order_parameter_sqrt3 rejects it (the docstring claim)
        latdef = obs_tri_lattice(4, 2)
        latdef.components[1][1] = true
        @test order_parameter_p2x2(latdef) ≈ sqrt(3) / 16 atol = 1e-12
        @test_throws ArgumentError order_parameter_sqrt3(latdef)

        # A short seeded canonical NS run records the Binder moments through
        # the moment-callback pattern from the docstring
        Random.seed!(4166)
        tri_ham1 = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
        walkers = [LatticeWalker(deepcopy(obs_tri_lattice(6, 2)),
                                 energy=0.0u"eV", iter=0) for _ in 1:20]
        # Fixed-N ladder at N = 6: place 6 particles per walker
        for w in walkers
            occ = vcat(fill(true, 6), fill(false, 18))
            shuffle!(occ)
            w.configuration.components[1] .= occ
        end
        ls = LatticeGasWalkers(walkers, tri_ham1; perturb_energy=1e-9)
        params = LatticeNestedSamplingParameters(mc_steps=30,
            energy_perturbation=1e-9, allowed_fail_count=100000)

        df, _, _ = nested_sampling(ls, params, Int64(200),
            MCRandomWalkClone(), obs_save;
            observables=[:psi => order_parameter_p2x2,
                         :psi2 => (cfg -> order_parameter_p2x2(cfg)^2),
                         :psi4 => (cfg -> order_parameter_p2x2(cfg)^4)])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "psi", "psi2", "psi4"]
        @test nrow(df) > 0
        @test all(0.0 .<= df.psi .<= 1 / 2)
        @test all(isapprox.(df.psi2, df.psi .^ 2; rtol=1e-12))
        @test all(isapprox.(df.psi4, df.psi .^ 4; rtol=1e-12))
    end

    @testset "observable hook: validation" begin
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:8]
        ls = LatticeGasWalkers(walkers, obs_ham)
        params = LatticeNestedSamplingParameters(mc_steps=10, energy_perturbation=1e-9)

        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=[:a => order_parameter_c2x2, :a => order_parameter_c2x2])
        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=[:emax => order_parameter_c2x2])
        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=Pair{Symbol,Function}[])
        @test_throws ArgumentError nested_sampling(ls, params, Int64(2),
            MCRandomWalkClone(), obs_save;
            observables=[:bad => (cfg -> "not a number")])
        # Non-finite probe values warn but do not throw
        @test_logs (:warn, r"non-finite") min_level=Base.CoreLogging.Warn match_mode=:any nested_sampling(
            ls, params, Int64(2), MCRandomWalkClone(), obs_save;
            observables=[:nanobs => (cfg -> NaN)])
        obs_cleanup()
    end

    # ================================================================
    # The pairing-guard test needs a step method that violates the cull
    # invariant; structs and method extensions must live at module top level
    @eval begin
        struct BrokenObsRoutine <: MCRoutine end
        function FreeBird.SamplingSchemes.nested_sampling_step!(
                liveset::LatticeGasWalkers,
                ns_params::NestedSamplingParameters,
                mc_routine::BrokenObsRoutine;
                ns_iteration::Int=0)
            # Claims an accepted iteration whose emax does not belong to the
            # pre-sorted worst walker
            return 1, liveset.walkers[1].energy + 1.0u"eV", liveset, ns_params
        end
    end

    @testset "observable hook: parallel rejection and pairing guard" begin
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                   for _ in 1:6]
        ls = LatticeGasWalkers(walkers, obs_ham)
        params = LatticeNestedSamplingParameters(mc_steps=5,
            energy_perturbation=1e-9)

        # Parallel/multi-cull routines cull several walkers per recorded row,
        # so they are rejected up front rather than corrupting the ledger
        @test_throws ArgumentError nested_sampling(ls, params, Int64(1),
            MCRandomWalkMaxEParallel(), obs_save;
            observables=[:psi => order_parameter_c2x2])

        # A step that reports a different walker than the pre-sorted worst
        # trips the bit-exact pairing guard
        @test_throws ErrorException nested_sampling(ls, params, Int64(1),
            BrokenObsRoutine(), obs_save;
            observables=[:psi => order_parameter_c2x2])
        obs_cleanup()
    end

    # ================================================================
    @testset "observable hook: exact dead-point pairing (igref route)" begin
        Random.seed!(42)
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:30]
        ls = LatticeGasWalkers(walkers, obs_ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(mc_steps=30, reference_fugacity=1.0)
        mc = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)

        df, final_ls, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(400), mc, obs_save;
            observables=[:n_check => count_N, :psi => order_parameter_c2x2])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "num_particles", "n_check", "psi"]
        @test eltype(df.n_check) == Float64
        @test eltype(df.psi) == Float64
        # The load-bearing check: the observable column, evaluated on the
        # captured culled walker, reproduces the independently recorded
        # num_particles column row by row over a full run
        @test df.n_check == Float64.(df.num_particles)
        @test all(0.0 .<= df.psi .<= 0.5)

        # Schema is unchanged when observables are not requested
        Random.seed!(42)
        walkers2 = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:30]
        ls2 = LatticeGasWalkers(walkers2, obs_ham; assign_energy=false)
        params2 = IdealGasReferencedGCNSParameters(mc_steps=30, reference_fugacity=1.0)
        df2, _, _ = ideal_gas_referenced_nested_sampling(
            ls2, params2, Int64(50), mc, obs_save)
        obs_cleanup()
        @test names(df2) == ["iter", "emax", "num_particles"]
    end

    # ================================================================
    @testset "observable hook: exact dead-point pairing (omega-sorted route)" begin
        Random.seed!(11)
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:20]
        ls = LatticeGasWalkers(walkers, obs_ham; assign_energy=false)
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=30, chemical_potential=-0.02)
        mc = MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3)

        df, _, _ = grand_canonical_nested_sampling(
            ls, gc_params, Int64(200), mc, obs_save;
            observables=[:n_check => count_N])
        obs_cleanup()

        @test names(df) == ["iter", "omega", "energy", "num_particles", "n_check"]
        @test df.n_check == Float64.(df.num_particles)
    end

    # ================================================================
    @testset "observable hook: exact dead-point pairing (canonical route)" begin
        Random.seed!(7)
        lat = obs_square_lattice(4, 4)
        walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0) for _ in 1:20]
        # Fixed-N ladder at N = 5: place 5 particles per walker
        for w in walkers
            occ = vcat(fill(true, 5), fill(false, 11))
            shuffle!(occ)
            w.configuration.components[1] .= occ
        end
        ls = LatticeGasWalkers(walkers, obs_ham; perturb_energy=1e-9)
        params = LatticeNestedSamplingParameters(mc_steps=30,
            energy_perturbation=1e-9, allowed_fail_count=100000)

        df, _, _ = nested_sampling(ls, params, Int64(200),
            MCRandomWalkClone(), obs_save;
            observables=[:e_check => (cfg -> interacting_energy(cfg, obs_ham).val),
                         :psi => order_parameter_c2x2])
        obs_cleanup()

        @test names(df) == ["iter", "emax", "e_check", "psi"]
        # N is conserved on the canonical route, so pair the recomputed
        # (unperturbed) energy against the recorded emax instead; they differ
        # only by the tie-breaking perturbation
        @test all(isapprox.(df.e_check, df.emax; atol=1e-8))
        @test all(0.0 .<= df.psi .<= 0.5)
    end

    # ================================================================
    @testset "igref stats: var_U, cov_UN, and observable averages" begin
        kb = 8.617333262e-5
        df = DataFrame(iter=[1, 2, 3], emax=[0.5, 0.3, 0.1],
                       num_particles=[2, 1, 0], a=[1.0, 2.0, 4.0])
        live_E = [0.05, 0.02]
        live_N = [1, 2]
        live_a = [8.0, 16.0]
        K, z0, M = 4, 1.0, 16
        ω0 = (K + 1) / K
        μs = [-0.02, 0.03]
        Ts = [300.0, 600.0]

        stats = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K;
            ω0=ω0, live_emax=live_E, live_numbers=live_N,
            observable_cols=[:a], live_observables=Dict(:a => live_a))

        # Independent linear-space reference for the closed-form weighted sums
        w0 = ω0 * (1 / (K + 1)) .* (K / (K + 1)) .^ df.iter
        wt = fill((K / (K + 1))^3 / K, 2)      # tail carries no ω0 factor
        w_all = vcat(w0, wt)
        E_all = vcat(df.emax, live_E)
        N_all = vcat(Float64.(df.num_particles), Float64.(live_N))
        a_all = vcat(df.a, live_a)
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            β = 1 / (kb * T)
            w = w_all .* exp.(β * μ .* N_all .- β .* E_all)   # z0 = 1
            sw = sum(w)
            u = sum(w .* E_all) / sw
            n = sum(w .* N_all) / sw
            @test stats.mean_U[i, j] ≈ u rtol = 1e-12
            @test stats.var_U[i, j] ≈ sum(w .* E_all .^ 2) / sw - u^2 rtol = 1e-10
            @test stats.cov_UN[i, j] ≈ sum(w .* E_all .* N_all) / sw - u * n rtol = 1e-10
            @test stats.observables[:a][i, j] ≈ sum(w .* a_all) / sw rtol = 1e-12
        end

        # Self-consistency: averaging the ledger's own columns reproduces
        # mean_N and mean_U
        # Mixed value types (Vector{Int} + Vector{Float64} typejoin to
        # Dict{Symbol,Vector}) pin the loosened Dict annotation: this is
        # exactly the documented self-consistency call
        stats2 = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K;
            ω0=ω0, live_emax=live_E, live_numbers=live_N,
            observable_cols=[:num_particles, :emax],
            live_observables=Dict(:num_particles => live_N, :emax => live_E))
        @test stats2.observables[:num_particles] ≈ stats2.mean_N rtol = 1e-12
        @test stats2.observables[:emax] ≈ stats2.mean_U rtol = 1e-12

        # Live-tail-only analysis: a zero-row ledger need not carry the
        # observable columns (mirrors the fixed-N empty-sector convention)
        df_empty = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[])
        st_tail = gc_thermodynamic_stats_ideal_ref(df_empty, M, z0, μs, Ts, K;
            ω0=ω0, live_emax=live_E, live_numbers=live_N,
            observable_cols=[:a], live_observables=Dict(:a => live_a))
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            β = 1 / (kb * T)
            w = exp.(β * μ .* Float64.(live_N) .- β .* live_E)
            @test st_tail.observables[:a][i, j] ≈
                  sum(w .* live_a) / sum(w) rtol = 1e-12
        end

        # Backward compatibility: pre-existing fields keep their order, new
        # fields are appended, positional destructuring still works
        stats3 = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K)
        @test isempty(stats3.observables)
        @test keys(stats3) == (:logXi, :mean_N, :var_N, :mean_U, :N_eff,
                               :var_U, :cov_UN, :observables)
        lX, mN, vN, mU, Ne = stats3
        @test lX == stats3.logXi && Ne == stats3.N_eff

        # Validation traps
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:nope])
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a, :a],
            live_emax=live_E, live_numbers=live_N,
            live_observables=Dict(:a => live_a))
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a],
            live_emax=live_E, live_numbers=live_N)
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a],
            live_emax=live_E, live_numbers=live_N,
            live_observables=Dict(:b => live_a))
        @test_throws DimensionMismatch gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K; observable_cols=[:a],
            live_emax=live_E, live_numbers=live_N,
            live_observables=Dict(:a => [1.0]))
        @test_throws ArgumentError gc_thermodynamic_stats_ideal_ref(
            df, M, z0, μs, Ts, K;
            live_observables=Dict(:a => live_a))
    end

    # ================================================================
    @testset "fixed-N stats: var_U, cov_UN, and observable averages" begin
        kb = 8.617333262e-5
        M = 4
        K = 3
        ω0 = (K + 1) / K
        N_values = [0, 1, 2]
        T_grid = [300.0, 600.0] .* u"K"
        μ_grid = [-0.05, 0.02] .* u"eV"

        df0 = DataFrame(iter=Int[], emax=Float64[])
        df1 = DataFrame(iter=[1, 2], emax=[0.3, 0.1], psi=[0.2, 0.4])
        df2 = DataFrame(iter=Int[], emax=Float64[])   # single-config convention
        live_E = [Float64[], [0.05, 0.01], [0.02, 0.02, 0.02]]
        live_psi = [Dict(:psi => [0.0]),
                    Dict(:psi => [0.6, 0.8]),
                    Dict(:psi => [0.5, 0.5, 0.5])]

        stats = gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, ω0=ω0, live_emax=live_E,
            observable_cols=[:psi], live_observables=live_psi)

        # Independent linear-space reference. Sector 1: two dead points plus
        # a two-entry tail (n_iters = 2). Sector 2: empty ladder, tail mass
        # exactly 1 split over three copies of the sector energy.
        w1 = vcat(ω0 * (1 / (K + 1)) .* (K / (K + 1)) .^ df1.iter,
                  fill(ω0 * (K / (K + 1))^2 / 2, 2))
        E1 = vcat(df1.emax, live_E[2])
        a1 = vcat(df1.psi, live_psi[2][:psi])
        binom = [1.0, 4.0, 6.0]                       # C(4, N)
        for (j, T) in enumerate([300.0, 600.0]), (k, μ) in enumerate([-0.05, 0.02])
            β = 1 / (kb * T)
            b1 = w1 .* exp.(-β .* E1)
            z1 = sum(b1)
            e1 = sum(b1 .* E1) / z1
            e1sq = sum(b1 .* E1 .^ 2) / z1
            p1 = sum(b1 .* a1) / z1
            z2 = exp(-β * 0.02)
            e2, e2sq, p2 = 0.02, 0.02^2, 0.5
            wN = binom .* [1.0, z1, z2] .* exp.(β * μ .* [0.0, 1.0, 2.0])
            sw = sum(wN)
            U = (wN[2] * e1 + wN[3] * e2) / sw
            Nbar = (wN[2] * 1 + wN[3] * 2) / sw
            @test stats.logXi[k, j] ≈ log(sw) rtol = 1e-10
            @test stats.mean_N[k, j] ≈ Nbar rtol = 1e-10
            @test stats.mean_U[k, j] ≈ U rtol = 1e-10
            @test stats.var_U[k, j] ≈
                  (wN[2] * e1sq + wN[3] * e2sq) / sw - U^2 rtol = 1e-8
            @test stats.cov_UN[k, j] ≈
                  (wN[2] * 1 * e1 + wN[3] * 2 * e2) / sw - U * Nbar rtol = 1e-8
            # The N = 0 sector contributes psi = 0 with weight wN[1]
            @test stats.observables[:psi][k, j] ≈
                  (wN[2] * p1 + wN[3] * p2) / sw rtol = 1e-10
        end

        # The N = 0 live_observables entry is USED (its live_emax counterpart
        # is ignored): with a nonzero Int-typed sentinel — which also pins the
        # loosened heterogeneous-Dict annotation — the empty sector's
        # contribution wN[1] * 10 must appear in <b>, and it dominates at the
        # strongly negative mu grid point
        df1b = DataFrame(iter=[1, 2], emax=[0.3, 0.1], psi=[0.2, 0.4],
                         b=[1.0, 3.0])
        live_b = [Dict{Symbol,Vector}(:psi => [0.0], :b => [10]),
                  Dict{Symbol,Vector}(:psi => [0.6, 0.8], :b => [5.0, 7.0]),
                  Dict{Symbol,Vector}(:psi => [0.5, 0.5, 0.5],
                                      :b => [2.0, 2.0, 2.0])]
        stats_b = gc_thermodynamic_stats_fixed_N(
            [df0, df1b, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, ω0=ω0, live_emax=live_E,
            observable_cols=[:psi, :b], live_observables=live_b)
        @test stats_b.observables[:psi] ≈ stats.observables[:psi] rtol = 1e-12
        for (j, T) in enumerate([300.0, 600.0]), (k, μ) in enumerate([-0.05, 0.02])
            β = 1 / (kb * T)
            b1v = w1 .* exp.(-β .* E1)
            z1 = sum(b1v)
            bb1 = sum(b1v .* vcat(df1b.b, live_b[2][:b])) / z1
            z2 = exp(-β * 0.02)
            wN = binom .* [1.0, z1, z2] .* exp.(β * μ .* [0.0, 1.0, 2.0])
            sw = sum(wN)
            @test stats_b.observables[:b][k, j] ≈
                  (wN[1] * 10.0 + wN[2] * bb1 + wN[3] * 2.0) / sw rtol = 1e-10
        end

        # Backward compatibility: field order preserved, new fields appended
        stats_nc = gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, ω0=ω0, live_emax=live_E)
        @test keys(stats_nc) == (:logXi, :mean_N, :var_N, :mean_U,
                                 :log_Z_N, :N_values, :var_U, :cov_UN,
                                 :observables)
        @test isempty(stats_nc.observables)
        @test stats_nc.logXi ≈ stats.logXi rtol = 1e-14

        # Validation traps
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi])
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, observable_cols=[:psi], live_observables=live_psi)
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=live_psi[1:2])
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=[Dict(:other => [0.0]), live_psi[2], live_psi[3]])
        @test_throws DimensionMismatch gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=[live_psi[1], Dict(:psi => [0.6]), live_psi[3]])
        df1_nocol = DataFrame(iter=[1, 2], emax=[0.3, 0.1])
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1_nocol, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, observable_cols=[:psi],
            live_observables=live_psi)
        @test_throws ArgumentError gc_thermodynamic_stats_fixed_N(
            [df0, df1, df2], N_values, M, μ_grid, T_grid;
            n_walkers=K, live_emax=live_E, live_observables=live_psi)
    end

    # ================================================================
    @testset "end-to-end: exact reference and route consistency (4x4, NN repulsion)" begin
        # Statistical NS checks: seed the global RNG (the random_seed field on
        # the NS parameters is not consumed). The igref tolerances were sized
        # at or above 3 sigma of three-seed scatter at this configuration
        # (max observed: |dlnXi| 0.23, |dN| 0.17, |dU| 0.0011, |dvar_U| 1e-4,
        # |dcov_UN| 0.0024, |dpsi| 0.014, all N_eff >= 248); the fixed-N
        # tolerances follow the existing exact-enumeration testset pattern.
        # The prior must sit near the target ensemble: with repulsive
        # couplings, z0 = 1 centers the Bernoulli prior at N = M/2 and the
        # mu grid at ln z0 = 0, and the ladder-depth formula below is the
        # full-compression depth M ln(1 + 1/z0).
        Random.seed!(4400)
        kb = 8.617333262e-5

        M = 16
        J = 0.1                      # eV, repulsive nearest-neighbor coupling
        ham = GenericLatticeHamiltonian(0.0, [J, 0.0], u"eV")
        lattice_at(N) = begin
            lat = obs_square_lattice(4, 4)
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end
        psi_from_occ(occ) = abs(sum(occ[s] ?
            (iseven(((s - 1) % 4) + ((s - 1) ÷ 4)) ? 1 : -1) : 0
            for s in eachindex(occ))) / length(occ)

        # Exact per-configuration (E, N, psi) from fixed-N enumeration of all
        # 2^16 configurations
        exact_E = Dict{Int,Vector{Float64}}()
        exact_psi = Dict{Int,Vector{Float64}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(lattice_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            exact_psi[N] = [psi_from_occ(cfg[1]) for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
        end

        function exact_stats(μ_val, T_val)
            β = 1.0 / (kb * T_val)
            lt = Float64[]; Ns = Float64[]; Es = Float64[]; Ps = Float64[]
            for N in 0:M, i in eachindex(exact_E[N])
                push!(lt, N * β * μ_val - β * exact_E[N][i])
                push!(Ns, N); push!(Es, exact_E[N][i]); push!(Ps, exact_psi[N][i])
            end
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            u = sum(w .* Es) / sw; n = sum(w .* Ns) / sw
            (logXi=maximum(lt) + log(sw), mean_N=n,
             var_N=sum(w .* Ns .^ 2) / sw - n^2, mean_U=u,
             var_U=sum(w .* Es .^ 2) / sw - u^2,
             cov_UN=sum(w .* Es .* Ns) / sw - u * n,
             psi=sum(w .* Ps) / sw)
        end

        μs = [-0.010, 0.000, 0.010]
        Ts = [300.0, 450.0]

        # ---- igref route with psi recording ----
        z0 = 1.0
        K_ig = 150
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig
        walkers = [LatticeWalker(deepcopy(lattice_at(0)), energy=0.0u"eV", iter=0)
                   for _ in 1:K_ig]
        ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        df_ig, final_ig, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(n_ig), MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3),
            obs_save; observables=[:psi => order_parameter_c2x2])
        obs_cleanup()
        live_E_ig = [w.energy.val for w in final_ig.walkers]
        live_N_ig = [sum(w.configuration.components[1]) for w in final_ig.walkers]
        live_psi_ig = [order_parameter_c2x2(w.configuration) for w in final_ig.walkers]

        st_ig = gc_thermodynamic_stats_ideal_ref(df_ig, M, z0, μs, Ts, K_ig;
            ω0=(K_ig + 1) / K_ig, live_emax=live_E_ig, live_numbers=live_N_ig,
            observable_cols=[:psi], live_observables=Dict(:psi => live_psi_ig))

        n_gated = 0
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            st_ig.N_eff[i, j] >= 200 || continue
            n_gated += 1
            ex = exact_stats(μ, T)
            @test isapprox(st_ig.logXi[i, j], ex.logXi, atol=0.7)
            @test isapprox(st_ig.mean_N[i, j], ex.mean_N, atol=0.5)
            @test isapprox(st_ig.mean_U[i, j], ex.mean_U, atol=0.01)
            @test isapprox(st_ig.var_U[i, j], ex.var_U, rtol=0.4, atol=0.001)
            @test isapprox(st_ig.cov_UN[i, j], ex.cov_UN, rtol=0.4, atol=0.01)
            @test isapprox(st_ig.observables[:psi][i, j], ex.psi, atol=0.05)
        end
        # The grid was chosen inside the z0 window, so most points must pass
        # the Kish gate
        @test n_gated >= 4

        # ---- fixed-N route with psi recording ----
        K_fn = 100
        n_fn = 1200
        dfs = Vector{DataFrame}(undef, M + 1)
        live_E_fn = Vector{Vector{Float64}}(undef, M + 1)
        live_psi_fn = Vector{Dict{Symbol,Vector{Float64}}}(undef, M + 1)

        dfs[1] = DataFrame(iter=Int[], emax=Float64[])
        live_E_fn[1] = Float64[]
        live_psi_fn[1] = Dict(:psi => [0.0])          # psi of the empty lattice

        for N in 1:(M - 1)
            walkers = [
                begin
                    lat = lattice_at(N)
                    generate_random_new_lattice_sample!(lat)
                    LatticeWalker(lat)
                end for _ in 1:K_fn]
            liveset = LatticeGasWalkers(walkers, ham, perturb_energy=1e-9)
            ns_params = LatticeNestedSamplingParameters(
                mc_steps=60, energy_perturbation=1e-9,
                allowed_fail_count=100_000)
            df, final_ls, _ = nested_sampling(
                liveset, ns_params, n_fn, MCRandomWalkClone(), obs_save;
                observables=[:psi => order_parameter_c2x2])
            dfs[N + 1] = df
            live_E_fn[N + 1] = [ustrip(u"eV", w.energy) for w in final_ls.walkers]
            live_psi_fn[N + 1] = Dict(
                :psi => [order_parameter_c2x2(w.configuration) for w in final_ls.walkers])
        end
        obs_cleanup()

        # N = M: single configuration; empty ladder + live tail with K copies
        dfs[M + 1] = DataFrame(iter=Int[], emax=Float64[])
        live_E_fn[M + 1] = fill(only(exact_E[M]), K_fn)
        live_psi_fn[M + 1] = Dict(:psi => fill(0.0, K_fn))   # full lattice: psi = 0

        st_fn = gc_thermodynamic_stats_fixed_N(dfs, collect(0:M), M,
            μs .* u"eV", Ts .* u"K";
            n_walkers=K_fn, n_cull=1, ω0=(K_fn + 1) / K_fn,
            live_emax=live_E_fn,
            observable_cols=[:psi], live_observables=live_psi_fn)

        for (j, T) in enumerate(Ts), (k, μ) in enumerate(μs)
            ex = exact_stats(μ, T)
            @test isapprox(st_fn.logXi[k, j], ex.logXi, atol=0.3)
            @test isapprox(st_fn.mean_N[k, j], ex.mean_N, rtol=0.05, atol=0.3)
            @test isapprox(st_fn.mean_U[k, j], ex.mean_U, rtol=0.05, atol=0.03)
            @test isapprox(st_fn.var_U[k, j], ex.var_U, rtol=0.3, atol=0.015)
            @test isapprox(st_fn.cov_UN[k, j], ex.cov_UN, rtol=0.3, atol=0.08)
            @test isapprox(st_fn.observables[:psi][k, j], ex.psi, atol=0.04)
        end

        # ---- route consistency where the igref window is healthy ----
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            st_ig.N_eff[i, j] >= 200 || continue
            @test isapprox(st_ig.observables[:psi][i, j],
                           st_fn.observables[:psi][i, j], atol=0.07)
            @test isapprox(st_ig.var_U[i, j], st_fn.var_U[i, j],
                           rtol=0.5, atol=0.005)
        end
    end

    # ================================================================
    @testset "end-to-end: triangular (3x3), NN repulsion, igref + omega-sorted" begin
        # First triangular coverage of both grand-canonical routes and the
        # first omega-sorted observable ledger. Tolerances sized at or above
        # 3 sigma of three-seed scatter at this configuration (max observed
        # over Kish-gated points: |dlnXi| 0.61, |dN| 0.18, |dU| 0.0038,
        # |dvar_U| 4.7e-4, |dcov_UN| 0.0075, |dpsi| 0.013, all gated
        # N_eff >= 206; sigma(lnXi) ≈ √(D/K) ≈ 0.29 at the full-compression
        # depth D = 18·ln 2, so atol = 0.9 ≈ 3.1σ).
        Random.seed!(11)
        kb = 8.617333262e-5

        M = 18
        J = 0.1                      # eV, repulsive nearest-neighbor coupling
        tri_ham = GenericLatticeHamiltonian(0.0, [J], u"eV")
        tri_at(N) = begin
            lat = obs_tri_lattice(3, 3)
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end
        # Independent psi: geometric labels from the positions (inverting
        # [t1 t2]) and the sublattice sum evaluated with complex arithmetic,
        # sharing no code with order_parameter_sqrt3's integer identity
        tri_labels = let tpl = obs_tri_lattice(3, 3), a = 1.0
            map(1:size(tpl.positions, 1)) do s
                x, y = tpl.positions[s, 1], tpl.positions[s, 2]
                n = round(Int, 2 * y / (sqrt(3) * a))
                m = round(Int, x / a - n / 2)
                mod(m - n, 3)
            end
        end
        tri_psi(occ) = abs(sum(occ[s] ? cis(2π * tri_labels[s] / 3) : 0.0im
                               for s in eachindex(occ))) / length(occ)

        # Exact per-configuration (E, N, psi) from fixed-N enumeration of all
        # 2^18 configurations
        exact_E = Dict{Int,Vector{Float64}}()
        exact_psi = Dict{Int,Vector{Float64}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(tri_at(N), tri_ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            exact_psi[N] = [tri_psi(cfg[1]) for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
        end

        function exact_stats_tri(μ_val, T_val)
            β = 1.0 / (kb * T_val)
            lt = Float64[]; Ns = Float64[]; Es = Float64[]; Ps = Float64[]
            for N in 0:M, i in eachindex(exact_E[N])
                push!(lt, N * β * μ_val - β * exact_E[N][i])
                push!(Ns, N); push!(Es, exact_E[N][i]); push!(Ps, exact_psi[N][i])
            end
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            u = sum(w .* Es) / sw; n = sum(w .* Ns) / sw
            (logXi=maximum(lt) + log(sw), mean_N=n, mean_U=u,
             var_U=sum(w .* Es .^ 2) / sw - u^2,
             cov_UN=sum(w .* Es .* Ns) / sw - u * n,
             psi=sum(w .* Ps) / sw)
        end

        # μ = 0.05 eV probes the Kish-window edge: it fails the gate at
        # 300 K (N_eff ≈ 20-50 across seeds) and passes at 450 K in most
        # seeds — the trust-window behavior the gate exists to detect
        μs = [-0.010, 0.000, 0.010, 0.050]
        Ts = [300.0, 450.0]

        # ---- igref route with psi recording ----
        z0 = 1.0
        K_ig = 150
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig
        walkers = [LatticeWalker(deepcopy(tri_at(0)), energy=0.0u"eV", iter=0)
                   for _ in 1:K_ig]
        ls = LatticeGasWalkers(walkers, tri_ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        df_ig, final_ig, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(n_ig), MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3),
            obs_save; observables=[:psi => order_parameter_sqrt3])
        obs_cleanup()
        live_E_ig = [w.energy.val for w in final_ig.walkers]
        live_N_ig = [sum(w.configuration.components[1]) for w in final_ig.walkers]
        live_psi_ig = [order_parameter_sqrt3(w.configuration) for w in final_ig.walkers]

        st_ig = gc_thermodynamic_stats_ideal_ref(df_ig, M, z0, μs, Ts, K_ig;
            ω0=(K_ig + 1) / K_ig, live_emax=live_E_ig, live_numbers=live_N_ig,
            observable_cols=[:psi], live_observables=Dict(:psi => live_psi_ig))

        # Ledger integrity: on this cell every configuration obeys psi <= 1/3
        @test names(df_ig) == ["iter", "emax", "num_particles", "psi"]
        @test all(0.0 .<= df_ig.psi .<= 1 / 3)

        n_gated = 0
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            st_ig.N_eff[i, j] >= 200 || continue
            n_gated += 1
            ex = exact_stats_tri(μ, T)
            @test isapprox(st_ig.logXi[i, j], ex.logXi, atol=0.9)
            @test isapprox(st_ig.mean_N[i, j], ex.mean_N, atol=0.5)
            @test isapprox(st_ig.mean_U[i, j], ex.mean_U, atol=0.01)
            @test isapprox(st_ig.var_U[i, j], ex.var_U, rtol=0.4, atol=0.001)
            @test isapprox(st_ig.cov_UN[i, j], ex.cov_UN, rtol=0.4, atol=0.015)
            @test isapprox(st_ig.observables[:psi][i, j], ex.psi, atol=0.025)
        end
        # The inner window (μ = 0, ±0.01) sits inside the z0 trust region
        @test n_gated >= 4

        # ---- omega-sorted route at μ = 0.05 eV with psi recording ----
        Random.seed!(11)
        K_gc = 100
        μ_gc = 0.05
        walkers_gc = [LatticeWalker(deepcopy(tri_at(0)), energy=0.0u"eV", iter=0)
                      for _ in 1:K_gc]
        ls_gc = LatticeGasWalkers(walkers_gc, tri_ham; assign_energy=false)
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=100, chemical_potential=μ_gc,
            energy_perturbation=1e-9, init_occupation_p=0.3)
        df_gc, _, _ = grand_canonical_nested_sampling(
            ls_gc, gc_params, Int64(3000),
            MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25), obs_save;
            observables=[:psi => order_parameter_sqrt3])
        obs_cleanup()

        @test names(df_gc) == ["iter", "omega", "energy", "num_particles", "psi"]
        @test all(0.0 .<= df_gc.psi .<= 1 / 3)

        β300 = 1.0 / (kb * 300.0)
        mean_E_gc, _, mean_N_gc = gc_thermodynamic_stats(df_gc, [β300], K_gc, μ_gc)
        ex300 = exact_stats_tri(μ_gc, 300.0)
        @test isapprox(mean_E_gc[1], ex300.mean_U, rtol=0.3)
        @test isapprox(mean_N_gc[1], ex300.mean_N, rtol=0.3)

        # Caller-side <psi> from the documented Ω-weight construction (the
        # legacy stats function carries no observable machinery by design;
        # this records the recipe): w_i ∝ Γ_i·exp(-β·Ω_i)
        ψ_at(β) = begin
            w = ωᵢ(df_gc.iter, K_gc) .*
                exp.(-β .* (df_gc.omega .- minimum(df_gc.omega)))
            sum(w .* df_gc.psi) / sum(w)
        end
        @test isapprox(ψ_at(β300), ex300.psi, atol=0.025)

        # ---- route consistency at the shared (μ, T) where igref is trusted ----
        β450 = 1.0 / (kb * 450.0)
        if st_ig.N_eff[4, 2] >= 200
            @test isapprox(st_ig.observables[:psi][4, 2], ψ_at(β450), atol=0.05)
        end
    end

    # ================================================================
    @testset "end-to-end: J1-J2 stripe model (4x4), igref + omega-sorted" begin
        # First multi-shell, mixed-sign end-to-end on both grand-canonical
        # routes: J1 = -0.1 eV (attractive first shell) with J2 = +0.1 eV
        # (repulsive second shell) selects the (2x1) superantiferromagnetic
        # stripe at half filling [Binder & Landau, PRB 21, 1941 (1980)].
        # Both shells are image-faithful at L = 4 (diagonal circumference
        # bound 2*sqrt(2) < 4), so the stock lattice fixture applies.
        # Statistical NS checks: seed the global RNG (the random_seed field
        # on the NS parameters is not consumed); the seeds are planted after
        # the enumeration block, which consumes RNG. Tolerances sized at or
        # above 3x the maximum three-seed deviation at this configuration
        # (max observed over the six Kish-gated points: |dlnXi| 0.201,
        # |dN| 0.218, |dU| 0.0223, |dvarU| 0.00256, |dcovUN| 0.0143 against
        # the exactly-zero mu = 0 covariance, |dpsi| 0.0137, all gated
        # N_eff >= 304), sanity-checked against sigma(lnXi) ~ sqrt(D/K)
        # ~ 0.27 at the full-compression depth D = 16 ln 2. The var_U and
        # cov_UN atol floors carry the 3x cushion unconditionally (3 x
        # 0.00256 -> 0.008; 3 x 0.0143 -> 0.045), independent of which grid
        # point produced the maximum; the rtol terms track the 600 K scale.
        kb = 8.617333262e-5

        M = 16
        ham = GenericLatticeHamiltonian(0.0, [-0.1, 0.1], u"eV")
        lattice_at(N) = begin
            lat = obs_square_lattice(4, 4)
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end
        # Independent stripe psi from the occupations: the two axial
        # period-2 Bragg amplitudes evaluated as (-1)^x / (-1)^y phase
        # sums, sharing no code with bragg_amplitude
        xy_of(s) = ((s - 1) % 4, ((s - 1) ÷ 4) % 4)
        psi_from_occ(occ) = begin
            zx = sum(occ[s] ? (-1.0)^xy_of(s)[1] : 0.0 for s in eachindex(occ))
            zy = sum(occ[s] ? (-1.0)^xy_of(s)[2] : 0.0 for s in eachindex(occ))
            sqrt(zx^2 + zy^2) / length(occ)
        end

        # Exact per-configuration (E, N, psi) from fixed-N enumeration of
        # all 2^16 configurations
        exact_E = Dict{Int,Vector{Float64}}()
        exact_psi = Dict{Int,Vector{Float64}}()
        exact_occ = Dict{Int,Vector{Vector{Bool}}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(lattice_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            exact_psi[N] = [psi_from_occ(cfg[1]) for cfg in df_exact.config]
            exact_occ[N] = [Vector{Bool}(cfg[1]) for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
        end

        # Ground-state pins, in units of the 0.1 eV coupling, over the
        # whole filling range
        @test [minimum(exact_E[N]) for N in 0:M] ≈
              0.1 .* [0, 0, -1, -2, -4, -4, -5, -6, -8,
                      -6, -5, -4, -4, -2, -1, 0, 0] atol = 1e-9
        # Shell coordinations and couplings cancel, 4*(-0.1) + 4*(+0.1) = 0,
        # so the model is exactly particle-hole symmetric: the full fixed-N
        # spectra mirror about N = 8
        @test all(isapprox(sort(exact_E[N]), sort(exact_E[M - N]); atol=1e-12)
                  for N in 0:M)
        # The N = 8 minimum, -0.8 eV, is attained by exactly the four
        # single-orientation stripe states (two row, two column), each with
        # full stripe order and no c(2x2) order
        gs8 = findall(e -> e <= minimum(exact_E[8]) + 1e-9, exact_E[8])
        @test length(gs8) == 4
        stripe_states = Set(
            Vector{Bool}([pred(xy_of(s)...) for s in 1:M])
            for pred in [(x, y) -> x % 2 == 0, (x, y) -> x % 2 == 1,
                         (x, y) -> y % 2 == 0, (x, y) -> y % 2 == 1])
        @test Set(exact_occ[8][i] for i in gs8) == stripe_states
        for i in gs8
            lat = lattice_at(0)
            lat.components[1] .= exact_occ[8][i]
            @test order_parameter_stripe(lat) ≈ 1 / 2 atol = 1e-12
            @test order_parameter_c2x2(lat) == 0.0
        end

        function exact_stats(μ_val, T_val)
            β = 1.0 / (kb * T_val)
            lt = Float64[]; Ns = Float64[]; Es = Float64[]; Ps = Float64[]
            for N in 0:M, i in eachindex(exact_E[N])
                push!(lt, N * β * μ_val - β * exact_E[N][i])
                push!(Ns, N); push!(Es, exact_E[N][i]); push!(Ps, exact_psi[N][i])
            end
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            u = sum(w .* Es) / sw; n = sum(w .* Ns) / sw
            (logXi=maximum(lt) + log(sw), mean_N=n, mean_U=u,
             var_U=sum(w .* Es .^ 2) / sw - u^2,
             cov_UN=sum(w .* Es .* Ns) / sw - u * n,
             psi=sum(w .* Ps) / sw)
        end

        # ... and the symmetric chemical potential is mu* = 0: from the
        # exact grand sums, <N> = 8 exactly at mu = 0 at any temperature
        @test exact_stats(0.0, 300.0).mean_N ≈ 8.0 atol = 1e-9
        @test exact_stats(0.0, 600.0).mean_N ≈ 8.0 atol = 1e-9

        # The stripe lobe is stable for mu in (-0.1, +0.1) eV; the grid
        # sits inside it, centered on the symmetric chemical potential
        μs = [-0.050, 0.000, 0.050]
        Ts = [300.0, 600.0]

        # ---- igref route with stripe recording ----
        Random.seed!(4501)
        z0 = 1.0
        K_ig = 150
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig
        walkers = [LatticeWalker(deepcopy(lattice_at(0)), energy=0.0u"eV", iter=0)
                   for _ in 1:K_ig]
        ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        df_ig, final_ig, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(n_ig), MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3),
            obs_save; observables=[:stripe => (cfg -> order_parameter_stripe(cfg))])
        obs_cleanup()
        live_E_ig = [w.energy.val for w in final_ig.walkers]
        live_N_ig = [sum(w.configuration.components[1]) for w in final_ig.walkers]
        live_psi_ig = [order_parameter_stripe(w.configuration) for w in final_ig.walkers]

        st_ig = gc_thermodynamic_stats_ideal_ref(df_ig, M, z0, μs, Ts, K_ig;
            ω0=(K_ig + 1) / K_ig, live_emax=live_E_ig, live_numbers=live_N_ig,
            observable_cols=[:stripe], live_observables=Dict(:stripe => live_psi_ig))

        # Ledger integrity: 0 <= Psi <= 1/2 row by row
        @test names(df_ig) == ["iter", "emax", "num_particles", "stripe"]
        @test all(0.0 .<= df_ig.stripe .<= 1 / 2)

        n_gated = 0
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            st_ig.N_eff[i, j] >= 200 || continue
            n_gated += 1
            ex = exact_stats(μ, T)
            @test isapprox(st_ig.logXi[i, j], ex.logXi, atol=0.7)
            @test isapprox(st_ig.mean_N[i, j], ex.mean_N, atol=0.7)
            @test isapprox(st_ig.mean_U[i, j], ex.mean_U, atol=0.07)
            @test isapprox(st_ig.var_U[i, j], ex.var_U, rtol=0.4, atol=0.008)
            @test isapprox(st_ig.cov_UN[i, j], ex.cov_UN, rtol=0.55, atol=0.045)
            @test isapprox(st_ig.observables[:stripe][i, j], ex.psi, atol=0.05)
        end
        # z0 = 1 centers the reference on the mu grid: all six points gated
        # at N_eff >= 304 across the calibration seeds (floor 200); one
        # point of slack retained, per the square end-to-end precedent
        @test n_gated >= 5

        # ---- omega-sorted route at mu = 0 eV with stripe recording ----
        Random.seed!(4501)
        K_gc = 100
        μ_gc = 0.0
        walkers_gc = [LatticeWalker(deepcopy(lattice_at(0)), energy=0.0u"eV", iter=0)
                      for _ in 1:K_gc]
        ls_gc = LatticeGasWalkers(walkers_gc, ham; assign_energy=false)
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=100, chemical_potential=μ_gc,
            energy_perturbation=1e-9, init_occupation_p=0.3)
        df_gc, _, _ = grand_canonical_nested_sampling(
            ls_gc, gc_params, Int64(3000),
            MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25), obs_save;
            observables=[:stripe => (cfg -> order_parameter_stripe(cfg))])
        obs_cleanup()

        @test names(df_gc) == ["iter", "omega", "energy", "num_particles", "stripe"]
        @test all(0.0 .<= df_gc.stripe .<= 1 / 2)

        β300 = 1.0 / (kb * 300.0)
        β600 = 1.0 / (kb * 600.0)
        mean_E_gc, _, mean_N_gc = gc_thermodynamic_stats(
            df_gc, [β300, β600], K_gc, μ_gc)
        ex300 = exact_stats(μ_gc, 300.0)
        ex600 = exact_stats(μ_gc, 600.0)
        @test isapprox(mean_E_gc[1], ex300.mean_U, rtol=0.3)
        @test isapprox(mean_N_gc[1], ex300.mean_N, rtol=0.3)
        @test isapprox(mean_E_gc[2], ex600.mean_U, rtol=0.3)
        @test isapprox(mean_N_gc[2], ex600.mean_N, rtol=0.3)

        # Caller-side <psi> from the documented Ω-weight construction (the
        # legacy stats function carries no observable machinery by design;
        # this re-records the recipe): w_i ∝ Γ_i·exp(-β·Ω_i)
        ψ_at(β) = begin
            w = ωᵢ(df_gc.iter, K_gc) .*
                exp.(-β .* (df_gc.omega .- minimum(df_gc.omega)))
            sum(w .* df_gc.stripe) / sum(w)
        end
        @test isapprox(ψ_at(β300), ex300.psi, atol=0.025)
        @test isapprox(ψ_at(β600), ex600.psi, atol=0.03)

        # ---- route consistency at the shared mu = 0 where igref is trusted ----
        for (j, β) in enumerate([β300, β600])
            if st_ig.N_eff[2, j] >= 200
                @test isapprox(st_ig.observables[:stripe][2, j], ψ_at(β), atol=0.05)
            end
        end
    end

    # ================================================================
    @testset "end-to-end: J1-J3 image-multiplicity model (4x4), igref" begin
        # First full-stack trip of image-multiplicity neighbor lists through
        # exact_enumeration, the igref sampler, and the stats layer, on the
        # nearest-neighbor-attraction plus isotropic third-neighbor-repulsion
        # model [Landau & Binder, PRB 31, 5946 (1985)]: J1 = -3u and
        # J3 = +2u with u = 1/30 eV, so every configuration energy is an
        # integer multiple of u. T = 0 bond count: the axial period-4 double
        # stripe pays (-3|J1| + 2*J3)/2 per site — three occupied NN entries
        # at J1/2 each plus two multiplicity-2 third-shell entries at J3/2 —
        # i.e. -1/12 eV per site, -0.6667 eV total at N = 8, strictly above
        # -0.8 eV; p(2x2) blocks and diagonal period-4 stripes pay -|J1| per
        # occupied site (two occupied NN entries) with zero third-shell
        # cost, winning for J3 > |J1|/2. Enumeration remains the arbiter of
        # every pinned value below.
        # Statistical NS checks: seed the global RNG (the random_seed field
        # on the NS parameters is not consumed); the seed is planted after
        # the enumeration block, which consumes RNG. Tolerances sized at or
        # above 3x the maximum three-seed deviation at this configuration
        # (max observed over the gated points: |dlnXi| 0.520, |dN| 0.335,
        # |dU| 0.0123, |dpsi| 0.0169, all gated N_eff >= 330); the ladder is
        # deeper than the J1-J2 one — D = 16 ln(1 + 1/0.075) ~ 42.6 nats,
        # sigma(lnXi) ~ sqrt(D/K) ~ 0.53 — so the lnXi gate is
        # proportionally looser (atol 1.6 ~ 3 sigma).
        kb = 8.617333262e-5

        M = 16
        uJ = 1 / 30                  # eV; J1 = -3*uJ, J3 = +2*uJ
        ham = GenericLatticeHamiltonian(0.0, [-0.1, 0.0, 1 / 15], u"eV")
        j3_lattice() = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.5, 2.1],
            components=[[false for _ in 1:16]],
            adsorptions=:full,
            image_multiplicity=true)
        lattice_at(N) = begin
            lat = j3_lattice()
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end

        # Structural pins: per-site shell coordinations exactly [4, 4, 4] —
        # the third shell reaches 2 distinct partner sites, each through 2
        # periodic images at L = 4
        tpl = j3_lattice()
        @test all(length.(tpl.neighbors[s]) == [4, 4, 4] for s in 1:M)
        @test all(length(unique(tpl.neighbors[s][3])) == 2 for s in 1:M)

        # Independent diagonal quadrature from the occupations, evaluated
        # as i^phase class sums, sharing no code with bragg_amplitude
        xy_of(s) = ((s - 1) % 4, ((s - 1) ÷ 4) % 4)
        psi_from_occ(occ) = begin
            zp = sum(occ[s] ? im^mod(xy_of(s)[1] + xy_of(s)[2], 4) : 0.0im
                     for s in eachindex(occ))
            zm = sum(occ[s] ? im^mod(xy_of(s)[1] - xy_of(s)[2], 4) : 0.0im
                     for s in eachindex(occ))
            sqrt(abs(zp)^2 + abs(zm)^2) / length(occ)
        end
        # The recorded observable: the diagonal k = (pi/2, +-pi/2)
        # quadrature composed caller-side from bragg_amplitude (the hook's
        # arbitrary-callable contract)
        diag_quad(cfg) = sqrt(bragg_amplitude(cfg, 1, 1)^2 +
                              bragg_amplitude(cfg, 1, -1)^2)

        # Exact per-configuration (E, N, psi) from fixed-N enumeration of
        # all 2^16 configurations, through the multiplicity neighbor lists
        exact_E = Dict{Int,Vector{Float64}}()
        exact_psi = Dict{Int,Vector{Float64}}()
        exact_occ = Dict{Int,Vector{Vector{Bool}}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(lattice_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            exact_psi[N] = [psi_from_occ(cfg[1]) for cfg in df_exact.config]
            exact_occ[N] = [Vector{Bool}(cfg[1]) for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
        end

        # Every configuration energy is an integer multiple of u
        @test all(all(isapprox(e / uJ, round(e / uJ); atol=1e-9)
                      for e in exact_E[N]) for N in 0:M)
        # Ground-state pins, in units of u, over the whole filling range
        @test [minimum(exact_E[N]) for N in 0:M] ≈
              uJ .* [0, 0, -3, -6, -12, -12, -15, -18, -24,
                     -22, -23, -24, -28, -26, -27, -28, -32] atol = 1e-9
        # The N = 8 minimum, -0.8 eV (-24u), has degeneracy exactly 16 —
        # the k = (pi/2, +-pi/2) manifold of p(2x2) blocks and diagonal
        # period-4 stripes; every member carries the full diagonal
        # quadrature and exactly no axial period-4 stripe order
        gs8 = findall(e -> e <= minimum(exact_E[8]) + 1e-9, exact_E[8])
        @test length(gs8) == 16
        for i in gs8
            lat = j3_lattice()
            lat.components[1] .= exact_occ[8][i]
            @test diag_quad(lat) ≈ sqrt(2) / 4 atol = 1e-12
            @test order_parameter_stripe(lat, period=4) == 0.0
            @test psi_from_occ(exact_occ[8][i]) ≈ sqrt(2) / 4 atol = 1e-12
        end

        # Grand-canonical stability from the E_min(N) convex hull: the
        # N = 8 phase is stable for mu in (-0.1, -1/30) eV, centered at
        # -1/15 eV — the igref grid below sits inside the lobe
        emin = [minimum(exact_E[N]) for N in 0:M]
        @test maximum((emin[9] - emin[N + 1]) / (8 - N) for N in 0:7) ≈
              -0.1 atol = 1e-9
        @test minimum((emin[N + 1] - emin[9]) / (N - 8) for N in 9:16) ≈
              -1 / 30 atol = 1e-9

        function exact_stats(μ_val, T_val)
            β = 1.0 / (kb * T_val)
            lt = Float64[]; Ns = Float64[]; Es = Float64[]; Ps = Float64[]
            for N in 0:M, i in eachindex(exact_E[N])
                push!(lt, N * β * μ_val - β * exact_E[N][i])
                push!(Ns, N); push!(Es, exact_E[N][i]); push!(Ps, exact_psi[N][i])
            end
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            u = sum(w .* Es) / sw; n = sum(w .* Ns) / sw
            (logXi=maximum(lt) + log(sw), mean_N=n, mean_U=u,
             psi=sum(w .* Ps) / sw)
        end

        μs = [-0.08, -1 / 15, -0.05]
        Ts = [300.0]

        # ---- igref route with the composed :psi_diag recording ----
        Random.seed!(4601)
        z0 = 0.075                   # ln z0 ~ beta*mu at the lobe center at 300 K
        K_ig = 150
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig
        walkers = [LatticeWalker(deepcopy(lattice_at(0)), energy=0.0u"eV", iter=0)
                   for _ in 1:K_ig]
        ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        df_ig, final_ig, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(n_ig), MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3),
            obs_save; observables=[:psi_diag => diag_quad])
        obs_cleanup()
        live_E_ig = [w.energy.val for w in final_ig.walkers]
        live_N_ig = [sum(w.configuration.components[1]) for w in final_ig.walkers]
        live_psi_ig = [diag_quad(w.configuration) for w in final_ig.walkers]

        st_ig = gc_thermodynamic_stats_ideal_ref(df_ig, M, z0, μs, Ts, K_ig;
            ω0=(K_ig + 1) / K_ig, live_emax=live_E_ig, live_numbers=live_N_ig,
            observable_cols=[:psi_diag],
            live_observables=Dict(:psi_diag => live_psi_ig))

        # Ledger integrity: the class-count bound (n0 - n2)^2 + (n1 - n3)^2
        # <= 32 per quadrature component caps the diagonal quadrature at 1/2
        @test names(df_ig) == ["iter", "emax", "num_particles", "psi_diag"]
        @test all(0.0 .<= df_ig.psi_diag .<= 1 / 2 + 1e-12)

        n_gated = 0
        for (i, μ) in enumerate(μs)
            st_ig.N_eff[i, 1] >= 200 || continue
            n_gated += 1
            ex = exact_stats(μ, 300.0)
            @test isapprox(st_ig.logXi[i, 1], ex.logXi, atol=1.6)
            @test isapprox(st_ig.mean_N[i, 1], ex.mean_N, atol=1.1)
            @test isapprox(st_ig.mean_U[i, 1], ex.mean_U, atol=0.04)
            @test isapprox(st_ig.observables[:psi_diag][i, 1], ex.psi, atol=0.06)
        end
        # ln z0 matched to the lobe center: all three points gated at
        # N_eff >= 330 across the calibration seeds (floor 200); one point
        # of slack retained, per the square end-to-end precedent
        @test n_gated >= 2
    end

    # ================================================================
    @testset "layer helpers: site_layers, layer_coverage, occupancy_profile, layer_field" begin
        slab = obs_slab_lattice(2, 2, 3)

        # lattice_positions ordering: basis innermost, dimension 1 fastest,
        # dimension 3 outermost, so layers are contiguous blocks
        @test site_layers(slab) == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]

        # Hand coverages: layer 1 full, one site of layer 2, layer 3 empty
        slab.components[1] .= false
        slab.components[1][1:4] .= true
        slab.components[1][6] = true
        @test layer_coverage(slab, 1) == 1.0
        @test layer_coverage(slab, 2) == 0.25
        @test layer_coverage(slab, 3) == 0.0
        @test occupancy_profile(slab) == [1.0, 0.25, 0.0]

        # mean(occupancy_profile) == N/M holds exactly in Float64: each
        # n_k/4 is a dyadic rational, their sum is exact, and the single
        # rounding of /3 equals the single rounding of N/12
        Random.seed!(29)
        for _ in 1:20
            occ = rand(Bool, 12)
            slab.components[1] .= occ
            @test mean(occupancy_profile(slab)) == sum(occ) / 12
        end

        # Two-site planar basis: n_basis enters the layer block size (8/layer)
        slab2b = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)],
            supercell_dimensions=(2, 2, 2),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:16]],
            adsorptions=:full,
            image_multiplicity=true)   # d1 = d2 = 2 wraps shell 1; see obs_slab_lattice
        @test site_layers(slab2b) == vcat(fill(1, 8), fill(2, 8))
        slab2b.components[1] .= false
        slab2b.components[1][1:3] .= true
        @test layer_coverage(slab2b, 1) == 3 / 8
        @test occupancy_profile(slab2b) == [3 / 8, 0.0]

        # d3 = 1 degenerates to the total coverage (documented, allowed)
        flat = obs_square_lattice(4, 4)
        flat.components[1] .= false
        flat.components[1][1:5] .= true
        @test site_layers(flat) == fill(1, 16)
        @test layer_coverage(flat, 1) == 5 / 16
        @test occupancy_profile(flat) == [5 / 16]

        # layer_field broadcast pins; element type preserved
        vals = [-0.27, -0.03375, -0.01]
        @test layer_field(obs_slab_lattice(2, 2, 3), vals) ==
              vcat(fill(-0.27, 4), fill(-0.03375, 4), fill(-0.01, 4))
        fu = layer_field(obs_slab_lattice(2, 2, 3), vals .* u"eV")
        @test eltype(fu) == typeof(1.0u"eV")
        @test fu == vcat(fill(-0.27, 4), fill(-0.03375, 4), fill(-0.01, 4)) .* u"eV"
        @test layer_field(flat, [0.5]) == fill(0.5, 16)
        # ... and the unitful profile feeds SiteFieldLatticeHamiltonian
        # directly, as the docstring recipe promises
        @test SiteFieldLatticeHamiltonian(
            GenericLatticeHamiltonian(0.0, [0.0], u"eV"), fu) isa
              SiteFieldLatticeHamiltonian

        # Layer bounds: 1 <= layer <= d3
        @test_throws ArgumentError layer_coverage(obs_slab_lattice(2, 2, 3), 0)
        @test_throws ArgumentError layer_coverage(obs_slab_lattice(2, 2, 3), 4)
        @test_throws ArgumentError layer_coverage(flat, 2)
        # layer_field length must equal d3
        @test_throws ArgumentError layer_field(obs_slab_lattice(2, 2, 3), [-0.27, -0.01])
        @test_throws ArgumentError layer_field(flat, [0.5, 0.5])
        # Non-planar basis: the z-block layer index is ill-defined, so all
        # four helpers refuse
        tilted = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
            supercell_dimensions=(2, 2, 2),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:16]],
            adsorptions=:full,
            image_multiplicity=true)
        @test_throws ArgumentError site_layers(tilted)
        @test_throws ArgumentError layer_coverage(tilted, 1)
        @test_throws ArgumentError occupancy_profile(tilted)
        @test_throws ArgumentError layer_field(tilted, [0.0, 0.0])
    end

    # ================================================================
    @testset "end-to-end: inhomogeneous Langmuir slab (2,2,3), site field, igref + theta ledger" begin
        # J = 0 with the inverse-cube layer profile eps_k = B/k^3,
        # B = -0.27 eV (the layer_field docstring example): the grand
        # ensemble factorizes exactly, Xi = prod_k (1 + x_k)^4 with
        # x_k = e^{beta(mu - eps_k)}, and theta_k = x_k/(1 + x_k), an
        # independent closed form sharing no code with the field energy
        # path. First end-to-end trip of a three-dimensional supercell
        # through enumeration, sampling, and the stats layer.
        # Statistical NS checks: seed the global RNG (the random_seed field
        # on the NS parameters is not consumed); the seed is planted after
        # the enumeration block, which consumes RNG. Tolerances sized at or
        # above 3x the maximum three-seed deviation at this configuration
        # (seeds 4701/4702/4703; max observed over the 18 gated points:
        # |dlnXi| 0.324, |dN| 0.136, |dU| 0.0034, |dtheta| 0.019, every
        # grid point gated in every seed at N_eff >= 549), sanity-checked
        # against sigma(lnXi) ~ sqrt(D/K) ~ 0.30 at the full-compression
        # depth D = 12 log1p(1/0.5) = 13.2 nats.
        kb = 8.617333262e-5

        M = 12
        eps_layer = [-0.27, -0.03375, -0.01]     # eV, B/k^3 with B = -0.27
        slab_at(N) = begin
            lat = obs_slab_lattice(2, 2, 3)
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end
        base = GenericLatticeHamiltonian(0.0, [0.0], u"eV")   # J = 0
        field = layer_field(obs_slab_lattice(2, 2, 3), eps_layer .* u"eV")
        ham = SiteFieldLatticeHamiltonian(base, field)

        # Independent per-site profile from hardcoded index arithmetic,
        # sharing no code with site_layers/layer_field
        site_eps = [eps_layer[(s - 1) ÷ 4 + 1] for s in 1:M]
        @test field == site_eps .* u"eV"   # welds the two oracles at the entrance

        # Exact per-configuration (E, N, n_k) from fixed-N enumeration of
        # all 2^12 = 4096 configurations; every energy must be the masked
        # field sum to 1e-12 eV, the exactness bound applied exhaustively
        exact_E = Dict{Int,Vector{Float64}}()
        exact_nk = Dict{Int,Vector{Vector{Int}}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(slab_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            occs = [Vector{Bool}(cfg[1]) for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
            for (e, occ) in zip(exact_E[N], occs)
                @test isapprox(e, sum(site_eps[occ]; init=0.0); atol=1e-12)
            end
            exact_nk[N] = [[sum(occ[4k-3:4k]) for k in 1:3] for occ in occs]
        end

        # Closed forms from the layer factorization
        langmuir(μ_val, T_val) = begin
            β = 1 / (kb * T_val)
            x = exp.(β .* (μ_val .- eps_layer))
            θ = x ./ (1 .+ x)
            (logXi=4 * sum(log1p, x), θ=θ,
             mean_N=4 * sum(θ), mean_U=4 * sum(eps_layer .* θ))
        end

        # Grand sums from the enumeration, assembled in the log domain
        function enum_stats(μ_val, T_val)
            β = 1 / (kb * T_val)
            lt = Float64[]; Ns = Float64[]; Es = Float64[]
            th = [Float64[] for _ in 1:3]
            for N in 0:M, i in eachindex(exact_E[N])
                push!(lt, N * β * μ_val - β * exact_E[N][i])
                push!(Ns, N); push!(Es, exact_E[N][i])
                for k in 1:3
                    push!(th[k], exact_nk[N][i][k] / 4)
                end
            end
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            (logXi=maximum(lt) + log(sw),
             mean_N=sum(w .* Ns) / sw, mean_U=sum(w .* Es) / sw,
             θ=[sum(w .* th[k]) / sw for k in 1:3])
        end

        # mu grid inside the igref trust window of z0 = 0.5 at 300 K; the
        # layers discriminate strongly across it (theta1 ~ 1, theta2 and
        # theta3 partially filled)
        μs = [-0.03, -0.015, 0.0]
        Ts = [300.0, 600.0]

        # Deterministic weld: enumeration == closed form at every grid
        # point (pure summation, house deterministic tolerance; the
        # extensive mean_N and mean_U carry a factor-M headroom)
        for T in Ts, μ in μs
            ex = langmuir(μ, T); en = enum_stats(μ, T)
            @test isapprox(en.logXi, ex.logXi; atol=1e-10)
            @test isapprox(en.mean_N, ex.mean_N; atol=1e-9)
            @test isapprox(en.mean_U, ex.mean_U; atol=1e-9)
            @test all(isapprox.(en.θ, ex.θ; atol=1e-10))
        end

        # ---- igref route recording theta1..theta3 (the sampler-hook
        # check rides this same run: identical fixture and ledger) ----
        Random.seed!(4701)
        z0 = 0.5                     # ln z0 = -0.69 centers the mu grid at 300 K
        K_ig = 150
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig
        walkers = [LatticeWalker(deepcopy(slab_at(0)), energy=0.0u"eV", iter=0)
                   for _ in 1:K_ig]
        ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        df_ig, final_ig, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(n_ig), MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3),
            obs_save;
            observables=[Symbol(:theta, k) => (cfg -> layer_coverage(cfg, k))
                         for k in 1:3])
        obs_cleanup()

        # Ledger integrity: columns, range, and the quarter-integer
        # lattice of 4-site-layer coverages
        @test names(df_ig) == ["iter", "emax", "num_particles",
                               "theta1", "theta2", "theta3"]
        for col in (df_ig.theta1, df_ig.theta2, df_ig.theta3)
            @test all(0.0 .<= col .<= 1.0)
            @test all(isinteger.(4 .* col))
        end

        live_E = [w.energy.val for w in final_ig.walkers]
        live_N = [sum(w.configuration.components[1]) for w in final_ig.walkers]
        live_th = Dict(Symbol(:theta, k) =>
                       [layer_coverage(w.configuration, k) for w in final_ig.walkers]
                       for k in 1:3)
        st = gc_thermodynamic_stats_ideal_ref(df_ig, M, z0, μs, Ts, K_ig;
            ω0=(K_ig + 1) / K_ig, live_emax=live_E, live_numbers=live_N,
            observable_cols=[:theta1, :theta2, :theta3], live_observables=live_th)

        # var_U and cov_UN are exercised by the existing end-to-end
        # testsets; the new code under test here is the field energy and
        # the theta observables, so those are what the gate asserts
        n_gated = 0
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            st.N_eff[i, j] >= 200 || continue
            n_gated += 1
            ex = langmuir(μ, T)
            @test isapprox(st.logXi[i, j], ex.logXi, atol=1.0)
            @test isapprox(st.mean_N[i, j], ex.mean_N, atol=0.5)
            @test isapprox(st.mean_U[i, j], ex.mean_U, atol=0.011)
            for k in 1:3
                @test isapprox(st.observables[Symbol(:theta, k)][i, j],
                               ex.θ[k], atol=0.06)
            end
        end
        # ln z0 centered on the 300 K grid: all six points gated at
        # N_eff >= 549 across the calibration seeds (floor 200); one
        # point of slack retained, per the end-to-end precedents
        @test n_gated >= 5
    end

    # ================================================================
    @testset "end-to-end: contact-substrate slab (3,3,2), igref + fixed-N + omega-sorted" begin
        # First three-dimensional end-to-end on the grand-canonical
        # routes, and the first three-route closure (ideal-gas-referenced,
        # fixed-N stitching, omega-sorted) on a single model: a bilayer
        # lattice gas with nearest-neighbor attraction and a strictly
        # contact substrate, the short-range-substrate layering regime of
        # [Pandit, Schick & Wortis, PRB 26, 5112 (1982)] in the
        # simple-cubic lattice-gas setting of [Binder & Landau, PRB 37,
        # 1745 (1988)]. The isotropic cell merges the in-plane and axial
        # nearest neighbors into one shell at distance 1.0; the in-plane
        # circumference 3 > 2 x 1.1 keeps the default minimum-image lists
        # bulk-tiled faithful, so construction is silent (contrast the
        # (2,2,3) fixtures, which need image multiplicity).
        # The igref row sits at 1450 K, the temperature that centers the
        # reference fugacity z0 = e^{beta mu} ~ 0.074 on the inter-step
        # midpoint: the reference walk's insertion acceptance scales as
        # z0 while ceiling-blocked deletions freeze the constrained walk,
        # so the route only mixes in N for z0 in the ~0.05-1 window every
        # shipped igref regression uses. Matching that window at the
        # staircase (mu ~ -0.33 eV) fixes the sampled row's temperature;
        # a 600 K row (z0 ~ 1.9e-3) reproducibly collapses onto the
        # full-lattice tail, <N> high by ~8. The crisp staircase physics
        # is carried at 300 K by the fixed-N and omega-sorted routes
        # below, which have no fugacity constraint.
        # Statistical NS checks: seed the global RNG (the random_seed
        # field on the NS parameters is not consumed); the seed is
        # planted after the enumeration block, which consumes RNG, and
        # re-planted at the start of each route so the three legs are
        # stream-independent. Tolerances sized at or above 3x the
        # maximum three-seed deviation at this configuration (seeds
        # 4801/4802/4803; igref maxima over the nine Kish-gated grid
        # evaluations: |dlnXi| 1.10, |dN| 1.10, |dU| 0.38,
        # |dtheta| 0.084, all gated N_eff >= 2767 on a floor of 200;
        # fixed-N maxima |dlnXi| 0.45, |dN| 0.18, |dU| 0.062,
        # |dtheta| 0.017; omega-route maxima |dE| 0.17, |dN| 0.49,
        # |dtheta| 7e-4 on the 300 K row and 0.041 (theta1) / 0.017
        # (theta2) on the 600 K row, gated per column; route-consistency
        # maxima |dtheta| 0.086 igref vs fixed-N and 0.134 igref vs
        # omega), sanity-checked against sigma(lnXi) ~
        # sqrt(D/K) ~ 0.57 at the reference-ladder depth
        # D = 18 log1p(1/z0) ~ 48 nats (worst igref |dlnXi| = 1.9
        # sigma).
        kb = 8.617333262e-5

        M = 18
        slab_a() = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(3, 3, 2),
            periodicity=(true, true, false),
            cutoff_radii=[1.1],
            components=[[false for _ in 1:18]],
            adsorptions=:full)
        lattice_at(N) = begin
            lat = slab_a()
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end
        field = layer_field(slab_a(), [-0.15, 0.0] .* u"eV")
        ham = SiteFieldLatticeHamiltonian(
            GenericLatticeHamiltonian(0.0, [-0.1], u"eV"), field)

        # Independent per-site field from hardcoded index arithmetic,
        # sharing no code with site_layers/layer_field: contact binding
        # on the nine layer-1 sites only
        site_eps = [s <= 9 ? -0.15 : 0.0 for s in 1:M]
        @test field == site_eps .* u"eV"

        # Faithfulness pins: one merged shell of five distinct partners
        # (4 in-plane + 1 axial) on every site, no image-doubled entries
        tpl = slab_a()
        @test all(length.(tpl.neighbors[s]) == [5] for s in 1:M)
        @test all(length(unique(tpl.neighbors[s][1])) == 5 for s in 1:M)
        @test site_layers(tpl) == vcat(fill(1, 9), fill(2, 9))

        # Exact per-configuration (E, N, n_k) from fixed-N enumeration
        # of all 2^18 configurations
        exact_E = Dict{Int,Vector{Float64}}()
        exact_nk = Dict{Int,Vector{Vector{Int}}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(lattice_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            exact_nk[N] = [[sum(cfg[1][1:9]), sum(cfg[1][10:18])]
                           for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
        end

        # T = 0 ladder. In-plane bond maxima on the 3x3 torus:
        # b(n) = [0, 0, 1, 3, 4, 6, 9, 11, 14, 18] for n = 0..9 (a full
        # row of three wraps into a 3-cycle; the complement identity
        # b(9 - k) = 18 - 4k + b_holes(k) fills in the rest). For N <= 9
        # the contact field wins and layer-1 packing is optimal:
        # E_min = -0.15 N - 0.1 b(N); for N = 9 + m the full first layer
        # is kept and each layer-2 site adds its axial bond:
        # E_min = -3.15 - 0.1 m - 0.1 b(m). Every value is a multiple of
        # 0.05 eV = gcd(|J|, eps_s); enumeration is the arbiter.
        emin = [minimum(exact_E[N]) for N in 0:M]
        @test emin ≈ 0.05 .* [0, -3, -8, -15, -20, -27, -36, -43, -52, -63,
                              -65, -69, -75, -79, -85, -93, -99, -107, -117] atol = 1e-9
        # The N = 9 ground state is unique and is exactly the full first
        # layer
        gs9 = findall(e -> e <= minimum(exact_E[9]) + 1e-9, exact_E[9])
        @test length(gs9) == 1
        @test exact_nk[9][only(gs9)] == [9, 0]

        # T = 0 staircase from the convex hull of E_min(N) (vertices
        # N = 0, 9, 18, so the isotherm has exactly two steps):
        # mu_1 = 2J - eps_s, the per-site cost of the first layer (18
        # in-plane bonds over 9 sites plus the contact field), and
        # mu_2 = 3J, the per-site cost of the second layer (two in-plane
        # bonds plus the axial bond, the T = 0 bulk-coexistence value).
        # Only the first step is substrate-shifted: the structural
        # signature of a strictly contact substrate.
        μ1h = maximum((emin[10] - emin[N + 1]) / (9 - N) for N in 0:8)
        μ2h = minimum((emin[N + 1] - emin[10]) / (N - 9) for N in 10:18)
        @test μ1h ≈ -0.35 atol = 1e-9
        @test μ2h ≈ -0.30 atol = 1e-9
        @test μ1h < μ2h
        μ_mid = (μ1h + μ2h) / 2

        # Grand sums on flat arrays built once (2^18 entries; the
        # isotherm scan below evaluates 51 grid points)
        flat_E = Float64[]; flat_N = Float64[]
        flat_t1 = Float64[]; flat_t2 = Float64[]
        for N in 0:M, i in eachindex(exact_E[N])
            push!(flat_E, exact_E[N][i]); push!(flat_N, N)
            push!(flat_t1, exact_nk[N][i][1] / 9)
            push!(flat_t2, exact_nk[N][i][2] / 9)
        end
        function exact_stats(μ_val, T_val)
            β = 1.0 / (kb * T_val)
            lt = β .* (μ_val .* flat_N .- flat_E)
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            (logXi=maximum(lt) + log(sw),
             mean_N=sum(w .* flat_N) / sw, mean_U=sum(w .* flat_E) / sw,
             θ1=sum(w .* flat_t1) / sw, θ2=sum(w .* flat_t2) / sw)
        end

        # Exact-isotherm plateau gates at 300 K: at the inter-step
        # midpoint mu = -0.325 the three cheapest excitations off the
        # layer-1-full state all cost Delta Omega = 0.225 eV (layer-1
        # hole 0.15 + 0.4 - 0.325; layer-2 adparticle -0.1 + 0.325;
        # collective layer-2 condensation 9 x (0.325 - 0.30)), so
        # exp(-0.225/kT) ~ 1.7e-4 and the 0.99/0.01 gates hold with
        # ~30x margin (they stay valid up to ~479 K, where the
        # multiplicity-weighted theta2 crosses 0.01). At the midpoint
        # the model also has an exact particle-hole x layer-swap
        # symmetry (2 mu - 5J + eps_s = 0 there), so <N> = 9 and
        # theta1 + theta2 = 1 identically at every temperature
        ex_p = exact_stats(μ_mid, 300.0)
        @test ex_p.θ1 >= 0.99
        @test ex_p.θ2 <= 0.01
        @test ex_p.mean_N ≈ 9.0 atol = 0.01
        @test exact_stats(-0.40, 300.0).mean_N <= 0.01
        @test exact_stats(-0.25, 300.0).mean_N >= 17.99
        # Staircase shape on a mu scan: layer ordering everywhere and
        # <N> nondecreasing (d<N>/dmu = beta var N >= 0 exactly)
        scan = [exact_stats(μ, 300.0) for μ in -0.45:0.005:-0.20]
        @test all(s.θ1 >= s.θ2 - 1e-12 for s in scan)
        @test all(diff([s.mean_N for s in scan]) .>= -1e-9)

        # mu grid straddling both steps (the outer points sit on the
        # smooth isotherm's risers); the 1450 K row is set by the z0
        # mixing window (header comment), well above the in-plane
        # condensation scale (2D Ising mapping kTc = 2.269 |J|/4,
        # Tc ~ 658 K), so this leg gates the reweighting machinery
        # against the exact grand sums rather than a sharp staircase
        μs = [-0.36, μ_mid, -0.29]
        Ts = [1450.0]

        # ---- igref route with theta recording ----
        Random.seed!(4801)
        z0 = exp(μ_mid / (kb * 1450.0))
        K_ig = 150
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig
        walkers = [LatticeWalker(deepcopy(lattice_at(0)), energy=0.0u"eV", iter=0)
                   for _ in 1:K_ig]
        ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        df_ig, final_ig, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(n_ig), MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3),
            obs_save;
            observables=[Symbol(:theta, k) => (cfg -> layer_coverage(cfg, k))
                         for k in 1:2])
        obs_cleanup()

        # Ledger integrity: columns, range, and the ninth-integer
        # lattice of 9-site-layer coverages
        @test names(df_ig) == ["iter", "emax", "num_particles", "theta1", "theta2"]
        for col in (df_ig.theta1, df_ig.theta2)
            @test all(0.0 .<= col .<= 1.0)
            @test all(isinteger.(9 .* col))
        end

        live_E_ig = [w.energy.val for w in final_ig.walkers]
        live_N_ig = [sum(w.configuration.components[1]) for w in final_ig.walkers]
        live_th_ig = Dict(Symbol(:theta, k) =>
                          [layer_coverage(w.configuration, k) for w in final_ig.walkers]
                          for k in 1:2)
        st_ig = gc_thermodynamic_stats_ideal_ref(df_ig, M, z0, μs, Ts, K_ig;
            ω0=(K_ig + 1) / K_ig, live_emax=live_E_ig, live_numbers=live_N_ig,
            observable_cols=[:theta1, :theta2], live_observables=live_th_ig)

        n_gated = 0
        for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
            st_ig.N_eff[i, j] >= 200 || continue
            n_gated += 1
            ex = exact_stats(μ, T)
            @test isapprox(st_ig.logXi[i, j], ex.logXi, atol=3.3)
            @test isapprox(st_ig.mean_N[i, j], ex.mean_N, atol=3.4)
            @test isapprox(st_ig.mean_U[i, j], ex.mean_U, atol=1.2)
            @test isapprox(st_ig.observables[:theta1][i, j], ex.θ1, atol=0.26)
            @test isapprox(st_ig.observables[:theta2][i, j], ex.θ2, atol=0.26)
            @test st_ig.observables[:theta1][i, j] >=
                  st_ig.observables[:theta2][i, j] - 0.05
        end
        # z0 centered on the mu grid: all three points gated at
        # N_eff >= 2767 across the calibration seeds (floor 200); one
        # point of slack retained, per the end-to-end precedents
        @test n_gated >= 2

        # ---- fixed-N route with theta recording, stitched over all 19
        # sectors ----
        Random.seed!(4801)
        K_fn = 100
        n_fn = 1300     # max sector depth ln C(18,9) ~ 10.8 nats -> 1.15 K x 10.8 ~ 1241
        Ts_fn = [300.0, 1450.0]
        dfs = Vector{DataFrame}(undef, M + 1)
        live_E_fn = Vector{Vector{Float64}}(undef, M + 1)
        live_th_fn = Vector{Dict{Symbol,Vector{Float64}}}(undef, M + 1)

        dfs[1] = DataFrame(iter=Int[], emax=Float64[])
        live_E_fn[1] = Float64[]
        live_th_fn[1] = Dict(:theta1 => [0.0], :theta2 => [0.0])   # empty lattice

        for N in 1:(M - 1)
            walkers = [
                begin
                    lat = lattice_at(N)
                    generate_random_new_lattice_sample!(lat)
                    LatticeWalker(lat)
                end for _ in 1:K_fn]
            liveset = LatticeGasWalkers(walkers, ham, perturb_energy=1e-9)
            ns_params = LatticeNestedSamplingParameters(
                mc_steps=60, energy_perturbation=1e-9,
                allowed_fail_count=100_000)
            df, final_ls, _ = nested_sampling(
                liveset, ns_params, n_fn, MCRandomWalkClone(), obs_save;
                observables=[Symbol(:theta, k) => (cfg -> layer_coverage(cfg, k))
                             for k in 1:2])
            dfs[N + 1] = df
            live_E_fn[N + 1] = [ustrip(u"eV", w.energy) for w in final_ls.walkers]
            live_th_fn[N + 1] = Dict(
                Symbol(:theta, k) => [layer_coverage(w.configuration, k)
                                      for w in final_ls.walkers] for k in 1:2)
        end
        obs_cleanup()

        # N = M: single configuration; empty ladder + live tail with K
        # copies; both layers full
        dfs[M + 1] = DataFrame(iter=Int[], emax=Float64[])
        live_E_fn[M + 1] = fill(only(exact_E[M]), K_fn)
        live_th_fn[M + 1] = Dict(:theta1 => fill(1.0, K_fn),
                                 :theta2 => fill(1.0, K_fn))

        @test names(dfs[10]) == ["iter", "emax", "theta1", "theta2"]

        st_fn = gc_thermodynamic_stats_fixed_N(dfs, collect(0:M), M,
            μs .* u"eV", Ts_fn .* u"K";
            n_walkers=K_fn, n_cull=1, ω0=(K_fn + 1) / K_fn,
            live_emax=live_E_fn,
            observable_cols=[:theta1, :theta2], live_observables=live_th_fn)

        # The 300 K row carries the crisp sampled staircase (exact
        # theta1 = 0.9997, theta2 = 3.4e-4 at the plateau midpoint), a
        # resolution outside the igref route's fugacity mixing window;
        # the fixed-N construction has no fugacity constraint. The
        # 1450 K row doubles as the igref route-consistency row.
        for (j, T) in enumerate(Ts_fn), (k, μ) in enumerate(μs)
            ex = exact_stats(μ, T)
            @test isapprox(st_fn.logXi[k, j], ex.logXi, atol=1.4)
            @test isapprox(st_fn.mean_N[k, j], ex.mean_N, rtol=0.05, atol=0.6)
            @test isapprox(st_fn.mean_U[k, j], ex.mean_U, rtol=0.05, atol=0.19)
            @test isapprox(st_fn.observables[:theta1][k, j], ex.θ1, atol=0.06)
            @test isapprox(st_fn.observables[:theta2][k, j], ex.θ2, atol=0.06)
        end

        # ---- omega-sorted route at the inter-step plateau ----
        Random.seed!(4801)
        K_gc = 100
        μ_gc = μ_mid
        walkers_gc = [LatticeWalker(deepcopy(lattice_at(0)), energy=0.0u"eV", iter=0)
                      for _ in 1:K_gc]
        ls_gc = LatticeGasWalkers(walkers_gc, ham; assign_energy=false)
        gc_params = GrandCanonicalNestedSamplingParameters(
            mc_steps=100, chemical_potential=μ_gc,
            energy_perturbation=1e-9, init_occupation_p=0.3)
        df_gc, _, _ = grand_canonical_nested_sampling(
            ls_gc, gc_params, Int64(3000),
            MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25), obs_save;
            observables=[Symbol(:theta, k) => (cfg -> layer_coverage(cfg, k))
                         for k in 1:2])
        obs_cleanup()

        @test names(df_gc) == ["iter", "omega", "energy", "num_particles",
                               "theta1", "theta2"]
        for col in (df_gc.theta1, df_gc.theta2)
            @test all(0.0 .<= col .<= 1.0)
            @test all(isinteger.(9 .* col))
        end

        β300 = 1.0 / (kb * 300.0)
        β600 = 1.0 / (kb * 600.0)
        β1450 = 1.0 / (kb * 1450.0)
        mean_E_gc, _, mean_N_gc = gc_thermodynamic_stats(
            df_gc, [β300, β600], K_gc, μ_gc)
        ex300 = exact_stats(μ_gc, 300.0)
        ex600 = exact_stats(μ_gc, 600.0)
        @test isapprox(mean_E_gc[1], ex300.mean_U, rtol=0.3)
        @test isapprox(mean_N_gc[1], ex300.mean_N, rtol=0.3)
        @test isapprox(mean_E_gc[2], ex600.mean_U, rtol=0.3)
        @test isapprox(mean_N_gc[2], ex600.mean_N, rtol=0.3)

        # Caller-side <theta_k> from the documented Ω-weight construction
        # (the legacy stats function carries no observable machinery by
        # design; this re-records the recipe): w_i ∝ Γ_i·exp(-β·Ω_i)
        θ_at(col, β) = begin
            w = ωᵢ(df_gc.iter, K_gc) .*
                exp.(-β .* (df_gc.omega .- minimum(df_gc.omega)))
            sum(w .* col) / sum(w)
        end
        @test isapprox(θ_at(df_gc.theta1, β300), ex300.θ1, atol=0.03)
        @test isapprox(θ_at(df_gc.theta2, β300), ex300.θ2, atol=0.02)
        @test isapprox(θ_at(df_gc.theta1, β600), ex600.θ1, atol=0.13)
        @test isapprox(θ_at(df_gc.theta2, β600), ex600.θ2, atol=0.06)
        @test θ_at(df_gc.theta1, β300) >= θ_at(df_gc.theta2, β300)
        @test θ_at(df_gc.theta1, β600) >= θ_at(df_gc.theta2, β600)

        # ---- route consistency ----
        # igref vs fixed-N on the shared 1450 K row, where the igref
        # window is healthy
        for (i, μ) in enumerate(μs)
            st_ig.N_eff[i, 1] >= 200 || continue
            @test isapprox(st_ig.observables[:theta1][i, 1],
                           st_fn.observables[:theta1][i, 2], atol=0.26)
            @test isapprox(st_ig.observables[:theta2][i, 1],
                           st_fn.observables[:theta2][i, 2], atol=0.26)
        end
        # igref vs omega-sorted at the shared plateau point, with the
        # omega ladder reweighted to the igref row's temperature
        if st_ig.N_eff[2, 1] >= 200
            @test isapprox(st_ig.observables[:theta1][2, 1],
                           θ_at(df_gc.theta1, β1450), atol=0.41)
            @test isapprox(st_ig.observables[:theta2][2, 1],
                           θ_at(df_gc.theta2, β1450), atol=0.41)
        end
    end

    # ================================================================
    @testset "end-to-end: tetragonal inverse-cube slab (2,2,3), igref + contact-biased channel" begin
        # The inverse-cube substrate staircase itself: nearest-neighbor
        # attraction [J_par, J_perp] on a tetragonal (2,2,3) slab with
        # the -eps_s/k^3 layer profile of [de Oliveira & Griffiths,
        # Surf. Sci. 71, 687 (1978)], whose T = 0 steps
        # mu_1 = 2 J_par - eps_s and mu_k = 2 J_par + J_perp - eps_s/k^3
        # are strictly monotone: the splitting of the upper steps that
        # part (a)'s strictly contact substrate collapses onto the
        # single bulk value 3J. First slab trip of image-multiplicity
        # neighbor lists and the interlayer-spacing geometry through
        # enumeration, sampling, and the stats layer, and the first
        # end-to-end exercise of the :contact biased insertion channel
        # and its per-move-type counters.
        # Statistical NS checks: seed the global RNG (the random_seed
        # field on the NS parameters is not consumed); the seed is
        # planted after the enumeration block, which consumes RNG, and
        # re-planted before the biased rerun. The rerun is identical in
        # parameters, not in stream: the biased channel consumes extra
        # RNG draws (the channel coin and the biased-site draw), so the
        # same seed diverges at the first insertion proposal; both runs
        # are gated against the same exact sums. The igref rows sit at
        # 1650 K, the temperature that centers z0 = e^{beta mu} ~ 0.073
        # on the second-plateau midpoint, inside the ~0.05-1 fugacity
        # window where the reference walk still mixes in N (see the
        # (3,3,2) testset's header); the staircase physics is
        # deterministic-only here. Tolerances sized at or above 3x the
        # maximum three-seed deviation over BOTH legs at this
        # configuration (seeds 4901/4902/4903; maxima over the eighteen
        # Kish-gated evaluations across both legs: |dlnXi| 0.68,
        # |dN| 0.88, |dU| 0.38, |dtheta| 0.13, all gated
        # N_eff >= 1623 on a floor of 200), sanity-checked against
        # sigma(lnXi) ~ sqrt(D/K) ~ 0.57 at the reference-ladder depth
        # D = 12 log1p(1/z0) ~ 32 nats (worst |dlnXi| = 1.2 sigma).
        kb = 8.617333262e-5

        M = 12
        eps_layer = [-0.27, -0.03375, -0.01]     # eV, -eps_s/k^3 with eps_s = 0.27
        slab_b() = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            interlayer_spacing=1.25,
            basis=[(0.0, 0.0, 0.0)],
            supercell_dimensions=(2, 2, 3),
            periodicity=(true, true, false),
            cutoff_radii=[1.1, 1.35],
            components=[[false for _ in 1:12]],
            adsorptions=:full,
            image_multiplicity=true)
        slab_at(N) = begin
            lat = slab_b()
            lat.components[1] .= false
            lat.components[1][1:N] .= true
            lat
        end
        base = GenericLatticeHamiltonian(0.0, [-0.1, -0.15], u"eV")   # [J_par, J_perp]
        field = layer_field(slab_b(), eps_layer .* u"eV")
        ham = SiteFieldLatticeHamiltonian(base, field)

        # Independent per-site profile from hardcoded index arithmetic,
        # sharing no code with site_layers/layer_field
        site_eps = [eps_layer[(s - 1) ÷ 4 + 1] for s in 1:M]
        @test field == site_eps .* u"eV"

        # Geometry pins. Shell 1 (in-plane, distance 1.0): the 2-wide
        # cell wraps, and multiplicity restores the bulk-tiled
        # coordination 4 as two distinct partners counted twice each.
        # Shell 2 (axial, distance 1.25): one bond on the surface layers,
        # two in the middle; the in-plane diagonal sqrt(2) and the axial
        # diagonal sqrt(1 + 1.25^2) = 1.60 both sit outside the cutoffs.
        tpl = slab_b()
        sl = site_layers(tpl)
        @test sl == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
        @test tpl.lattice_vectors[3, 3] == 1.25
        @test all(length(tpl.neighbors[s][1]) == 4 for s in 1:M)
        @test all(length(unique(tpl.neighbors[s][1])) == 2 for s in 1:M)
        @test all(length(tpl.neighbors[s][2]) == (sl[s] == 2 ? 2 : 1) for s in 1:M)

        # Exact per-configuration (E, N, n_k) from fixed-N enumeration
        # of all 2^12 = 4096 configurations
        exact_E = Dict{Int,Vector{Float64}}()
        exact_nk = Dict{Int,Vector{Vector{Int}}}()
        for N in 0:M
            df_exact, _ = exact_enumeration(slab_at(N), ham)
            exact_E[N] = [ustrip(u"eV", e) for e in df_exact.energy]
            exact_nk[N] = [[sum(cfg[1][4k-3:4k]) for k in 1:3]
                           for cfg in df_exact.config]
            @test length(exact_E[N]) == binomial(M, N)
        end

        # T = 0 ladder in the bulk-tiled doubled-bond convention (a full
        # 2x2 layer carries 8 in-plane entries = per-site 2 J_par, the
        # infinite-lattice value the multiplicity lists exist to
        # restore). In-plane entry counts b(n) on the wrapped 2x2 layer:
        # [0, 0, 2, 4, 8]. Layers fill bottom-up; each site in layer
        # k >= 2 adds one axial bond to the full layer below:
        #   N = 0..4:  E = -0.27 N - 0.1 b(N)
        #   N = 4+m:   E = E(4) - (0.03375 + 0.15) m - 0.1 b(m)
        #   N = 8+m:   E = E(8) - (0.01 + 0.15) m - 0.1 b(m)
        # No clean single unit exists (gcd granularity 1/800 eV), so the
        # decimal literals are pinned; enumeration is the arbiter.
        emin = [minimum(exact_E[N]) for N in 0:M]
        @test emin ≈ [0.0, -0.27, -0.74, -1.21, -1.88,
                      -2.06375, -2.4475, -2.83125, -3.415,
                      -3.575, -3.935, -4.295, -4.855] atol = 1e-9
        # The one- and two-layer film ground states are unique and are
        # exactly the bottom-up layer fillings
        gs4 = findall(e -> e <= minimum(exact_E[4]) + 1e-9, exact_E[4])
        @test length(gs4) == 1
        @test exact_nk[4][only(gs4)] == [4, 0, 0]
        gs8 = findall(e -> e <= minimum(exact_E[8]) + 1e-9, exact_E[8])
        @test length(gs8) == 1
        @test exact_nk[8][only(gs8)] == [4, 4, 0]

        # Three-step staircase from the hull (vertices N = 0, 4, 8, 12):
        # mu_1 = 2 J_par - eps_s; mu_k = 2 J_par + J_perp - eps_s/k^3
        # for k = 2, 3 -- strictly monotone, the inverse-cube splitting
        # this fixture exists to pin
        μ1h = maximum((emin[5] - emin[N + 1]) / (4 - N) for N in 0:3)
        μ2h = maximum((emin[9] - emin[N + 1]) / (8 - N) for N in 0:7)
        μ3h = minimum((emin[N + 1] - emin[9]) / (N - 8) for N in 9:12)
        @test μ1h ≈ -0.47 atol = 1e-9
        @test μ2h ≈ -0.38375 atol = 1e-9
        @test μ3h ≈ -0.36 atol = 1e-9
        @test μ1h < μ2h < μ3h

        # Grand sums on flat arrays built once
        flat_E = Float64[]; flat_N = Float64[]
        flat_th = [Float64[] for _ in 1:3]
        for N in 0:M, i in eachindex(exact_E[N])
            push!(flat_E, exact_E[N][i]); push!(flat_N, N)
            for k in 1:3
                push!(flat_th[k], exact_nk[N][i][k] / 4)
            end
        end
        function exact_stats(μ_val, T_val)
            β = 1.0 / (kb * T_val)
            lt = β .* (μ_val .* flat_N .- flat_E)
            w = exp.(lt .- maximum(lt)); sw = sum(w)
            (logXi=maximum(lt) + log(sw),
             mean_N=sum(w .* flat_N) / sw, mean_U=sum(w .* flat_E) / sw,
             θ=[sum(w .* flat_th[k]) / sw for k in 1:3])
        end

        # Exact-isotherm plateau gates. The resolving temperature is set
        # by the collective competitor, not the single-site one: on the
        # second plateau the nearest competing states are N = 4 and
        # N = 12 at Delta Omega = 4 x (0.02375/2) = 0.0475 eV
        # (single-site defects cost 0.2119 eV), so the 0.99/0.01 form
        # needs kT <= 0.0475/ln 100, i.e. T <= 120 K; at 300 K theta_3
        # reaches 0.121 on this plateau. Shipped at 100 K:
        # theta_2 = 0.9960, theta_3 = 0.0040, a 2.5x margin.
        T_iso = 100.0
        μ_c1 = (μ1h + μ2h) / 2     # first-plateau midpoint, -0.426875
        μ_c2 = (μ2h + μ3h) / 2     # second-plateau midpoint, -0.371875
        ex1 = exact_stats(μ_c1, T_iso)
        @test ex1.θ[1] >= 0.99
        @test ex1.θ[2] <= 0.01
        @test ex1.mean_N ≈ 4.0 atol = 0.01
        ex2 = exact_stats(μ_c2, T_iso)
        @test ex2.θ[1] >= 0.99
        @test ex2.θ[2] >= 0.99
        @test ex2.θ[3] <= 0.01
        @test ex2.mean_N ≈ 8.0 atol = 0.01
        @test exact_stats(-0.50, T_iso).mean_N <= 0.01
        @test exact_stats(-0.33, T_iso).mean_N >= 11.99
        scan = [exact_stats(μ, T_iso) for μ in -0.50:0.005:-0.315]
        @test all(s.θ[1] >= s.θ[2] - 1e-12 && s.θ[2] >= s.θ[3] - 1e-12
                  for s in scan)
        @test all(diff([s.mean_N for s in scan]) .>= -1e-9)

        # mu grid straddling the second and third steps. The 1650 K row
        # is set by the z0 mixing window (header comment); the sampled
        # gates are against the exact grand sums at the same (mu, T),
        # and the 0.27 eV contact field is ~2 kT there, so the layer
        # profile stays column-discriminating
        μs = [-0.40, μ_c2, -0.34]
        Ts = [1650.0]
        z0 = exp(μ_c2 / (kb * 1650.0))
        K_ig = 100
        n_ig = ceil(Int, 1.15 * K_ig * M * log1p(1 / z0)) + 2 * K_ig

        # One igref run shape for both legs: identical parameters by
        # construction, only the move routine differs
        run_igref(mc) = begin
            walkers = [LatticeWalker(deepcopy(slab_at(0)), energy=0.0u"eV", iter=0)
                       for _ in 1:K_ig]
            ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
            params = IdealGasReferencedGCNSParameters(
                mc_steps=100, reference_fugacity=z0, energy_perturbation=1e-9,
                allowed_fail_count=100_000)
            df, final_ls, params_out = ideal_gas_referenced_nested_sampling(
                ls, params, Int64(n_ig), mc, obs_save;
                observables=[Symbol(:theta, k) => (cfg -> layer_coverage(cfg, k))
                             for k in 1:3])
            obs_cleanup()
            (df, final_ls, params_out)
        end

        # One gate block for both legs: same ledger pins, same Kish
        # gate, same tolerances
        function gate_leg(df, final_ls)
            @test names(df) == ["iter", "emax", "num_particles",
                                "theta1", "theta2", "theta3"]
            for col in (df.theta1, df.theta2, df.theta3)
                @test all(0.0 .<= col .<= 1.0)
                @test all(isinteger.(4 .* col))
            end
            live_E = [w.energy.val for w in final_ls.walkers]
            live_N = [sum(w.configuration.components[1]) for w in final_ls.walkers]
            live_th = Dict(Symbol(:theta, k) =>
                           [layer_coverage(w.configuration, k) for w in final_ls.walkers]
                           for k in 1:3)
            st = gc_thermodynamic_stats_ideal_ref(df, M, z0, μs, Ts, K_ig;
                ω0=(K_ig + 1) / K_ig, live_emax=live_E, live_numbers=live_N,
                observable_cols=[:theta1, :theta2, :theta3],
                live_observables=live_th)
            n_gated = 0
            for (j, T) in enumerate(Ts), (i, μ) in enumerate(μs)
                st.N_eff[i, j] >= 200 || continue
                n_gated += 1
                ex = exact_stats(μ, T)
                @test isapprox(st.logXi[i, j], ex.logXi, atol=2.1)
                @test isapprox(st.mean_N[i, j], ex.mean_N, atol=2.7)
                @test isapprox(st.mean_U[i, j], ex.mean_U, atol=1.2)
                for k in 1:3
                    @test isapprox(st.observables[Symbol(:theta, k)][i, j],
                                   ex.θ[k], atol=0.4)
                end
                @test st.observables[:theta1][i, j] >=
                      st.observables[:theta2][i, j] - 0.05
                @test st.observables[:theta2][i, j] >=
                      st.observables[:theta3][i, j] - 0.05
            end
            # z0 centered on the mu grid: all three points gated at
            # N_eff >= 1623 across the calibration seeds and both legs
            # (floor 200); one point of slack retained, per the
            # end-to-end precedents
            @test n_gated >= 2
            st
        end

        # ---- igref route, uniform insertions ----
        Random.seed!(4901)
        df_u, final_u, params_u = run_igref(
            MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3))
        gate_leg(df_u, final_u)
        # The uniform run never touches the biased channel
        @test get(params_u.move_stats, :insert_biased_attempted, 0) == 0

        # ---- the identical run under the :contact biased channel ----
        # bias_shells = 2 spans both bond classes (in-plane and axial),
        # so a contact site is any empty site touching the film
        Random.seed!(4901)
        df_b, final_b, params_b = run_igref(MCGrandCanonicalMoves(
            p_move=0.4, p_insert=0.3, p_bias=0.5,
            bias_predicate=:contact, bias_shells=2))
        gate_leg(df_b, final_b)
        # First end-to-end trip of the biased channel: the run totals
        # must show the 0.5 split really split, and biased proposals
        # really accepted
        @test params_b.move_stats[:insert_biased_attempted] > 0
        @test params_b.move_stats[:insert_biased_accepted] > 0
        @test params_b.move_stats[:insert_uniform_attempted] > 0
        @test params_b.move_stats[:insert_biased_accepted] <=
              params_b.move_stats[:insert_biased_attempted]
    end

    # ================================================================
    @testset "end-to-end: triangular eight-shell + trio cluster expansion, igref" begin
        # Independent geometry for the triangular cells: explicit-image torus
        # distances from the stored positions, no FreeBird neighbor or
        # embedding code
        tce_sq3 = sqrt(3.0)
        function tce_tdist(p, q, Lx, Ly)
            best = Inf
            for ix in -1:1, iy in -1:1
                best = min(best, hypot(p[1] - q[1] + ix * Lx, p[2] - q[2] + iy * Ly))
            end
            return best
        end

        # Structural pins on the 18-site (3, 3, 1) cell: the NN-triangle
        # motif counts 42 embeddings under the torus convention (36 faces
        # plus 6 winding three-cycles) and the linear trio 36 (aliasing),
        # then the full 2^18 spectrum of the pair + two-trio expansion is
        # checked configuration-resolved against an independent evaluation
        tri18 = obs_tri_lattice(3, 3)
        tri18_tris = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            tri18, [1.0, 1.0, 1.0]; expected_count=42)
        tri18_lins = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            tri18, [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])
        @test length(tri18_lins) == 36
        tce_ham18 = ClusterLatticeHamiltonian(
            GenericLatticeHamiltonian(0.0, [0.1], u"eV"),
            [ClusterInteraction(0.05u"eV", tri18_tris),
             ClusterInteraction(-0.02u"eV", tri18_lins)])
        p18 = [(tri18.positions[s, 1], tri18.positions[s, 2]) for s in 1:18]
        pairs18 = [(i, j) for i in 1:18 for j in (i+1):18
                   if abs(tce_tdist(p18[i], p18[j], 3.0, 3 * tce_sq3) - 1.0) < 1e-9]
        @test length(pairs18) == 54
        probe18 = deepcopy(tri18)
        worst18 = 0.0
        for mask in 0:(2^18-1)
            occ = probe18.components[1]
            for s in 1:18
                occ[s] = ((mask >> (s - 1)) & 1) == 1
            end
            e_own = 0.1 * count(pr -> occ[pr[1]] && occ[pr[2]], pairs18) +
                    0.05 * count(t -> occ[t[1]] && occ[t[2]] && occ[t[3]], tri18_tris) -
                    0.02 * count(t -> occ[t[1]] && occ[t[2]] && occ[t[3]], tri18_lins)
            worst18 = max(worst18, abs(ustrip(u"eV", interacting_energy(probe18, tce_ham18)) - e_own))
        end
        @test worst18 <= 1e-12

        # End-to-end cell: (4, 2, 1) with image_multiplicity (the bulk-tiled
        # convention restores the coordination sequence (6, 6, 6, 12, 6, 6,
        # 12, 6) on the small cell), eight pair shells plus both trio motifs
        # enumerated on this cell with independently counted expected_count.
        # The NN-triangle enumeration is faithful here (shortest circumference
        # 2sqrt(3) > 3·1); the linear trio wraps and warns.
        tce_lat = MLattice{1,TriangularLattice}(
            lattice_constant=1.0,
            supercell_dimensions=(4, 2, 1),
            periodicity=(true, true, false),
            cutoff_radii=[1.05, 1.75, 2.05, 2.7, 3.05, 3.5, 3.65, 4.05],
            components=[[false for _ in 1:16]],
            adsorptions=:full,
            image_multiplicity=true)
        @test all(nb -> length.(nb) == [6, 6, 6, 12, 6, 6, 12, 6], tce_lat.neighbors)
        p16 = [(tce_lat.positions[s, 1], tce_lat.positions[s, 2]) for s in 1:16]
        function tce_trios(dists)
            found = NTuple{3,Int}[]
            for a in 1:16, b in (a+1):16, c in (b+1):16
                ds = sort([tce_tdist(p16[a], p16[b], 4.0, 2 * tce_sq3),
                           tce_tdist(p16[a], p16[c], 4.0, 2 * tce_sq3),
                           tce_tdist(p16[b], p16[c], 4.0, 2 * tce_sq3)])
                all(abs.(ds .- dists) .< 1e-9) && push!(found, (a, b, c))
            end
            return found
        end
        own_tris = tce_trios([1.0, 1.0, 1.0])
        own_lins = tce_trios([1.0, 1.0, 2.0])
        tce_tris = enumerate_motif_embeddings(tce_lat, [1.0, 1.0, 1.0];
            expected_count=length(own_tris))
        tce_lins = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            tce_lat, [1.0, 1.0, 2.0]; expected_count=length(own_lins))
        @test sort(collect(tce_tris)) == sort(own_tris)
        @test sort(collect(tce_lins)) == sort(own_lins)
        tce_J8 = [0.10, -0.02, 0.01, -0.005, 0.004, -0.003, 0.002, -0.001]
        tce_ham = ClusterLatticeHamiltonian(
            GenericLatticeHamiltonian(0.0, tce_J8, u"eV"),
            [ClusterInteraction(0.05u"eV", tce_tris),
             ClusterInteraction(-0.02u"eV", tce_lins)])

        # Independent multiplicity-weighted bond multiset from explicit image
        # enumeration, then the full 2^16 own-energy table and the
        # configuration-resolved spectrum check; the E_min(N) ladder's ends
        # pin the empty lattice at exactly zero and the full lattice at the
        # bulk-tiled closed form
        tce_radii = [1.05, 1.75, 2.05, 2.7, 3.05, 3.5, 3.65, 4.05]
        bonds16 = Tuple{Int,Int,Int}[]
        for i in 1:16, j in 1:16
            for ix in -3:3, iy in -3:3
                (j == i && ix == 0 && iy == 0) && continue
                d = hypot(p16[i][1] - p16[j][1] + ix * 4.0,
                          p16[i][2] - p16[j][2] + iy * 2 * tce_sq3)
                d > tce_radii[end] && continue
                sh = findfirst(r -> d <= r, tce_radii)
                sh === nothing || push!(bonds16, (i, j, sh))
            end
        end
        occ16 = fill(false, 16)
        e_own16 = zeros(2^16)
        for mask in 0:(2^16-1)
            for s in 1:16
                occ16[s] = ((mask >> (s - 1)) & 1) == 1
            end
            e = 0.0
            for (i, j, sh) in bonds16
                (occ16[i] && occ16[j]) && (e += tce_J8[sh] / 2)
            end
            e += 0.05 * count(t -> occ16[t[1]] && occ16[t[2]] && occ16[t[3]], tce_tris)
            e -= 0.02 * count(t -> occ16[t[1]] && occ16[t[2]] && occ16[t[3]], tce_lins)
            e_own16[mask+1] = e
        end
        probe16 = deepcopy(tce_lat)
        psi16 = zeros(2^16)
        n16 = zeros(Int, 2^16)
        worst16 = 0.0
        for mask in 0:(2^16-1)
            occ = probe16.components[1]
            for s in 1:16
                occ[s] = ((mask >> (s - 1)) & 1) == 1
            end
            worst16 = max(worst16, abs(ustrip(u"eV", interacting_energy(probe16, tce_ham)) - e_own16[mask+1]))
            psi16[mask+1] = order_parameter_p2x2(probe16)
            n16[mask+1] = count_ones(mask)
        end
        @test worst16 <= 1e-12
        emin16 = [minimum(e_own16[mask+1] for mask in 0:(2^16-1) if count_ones(mask) == N)
                  for N in 0:16]
        @test emin16[1] == 0.0
        tce_zbulk = [6, 6, 6, 12, 6, 6, 12, 6]
        tce_full = 16 / 2 * sum(tce_zbulk .* tce_J8) +
                   0.05 * length(tce_tris) - 0.02 * length(tce_lins)
        @test abs(emin16[17] - tce_full) < 1e-12

        # Constructed-configuration pins on the perfect p(2x2): axial
        # coordinates from the positions, occupied where both are even
        tce_ij = map(1:16) do s
            x, y = tce_lat.positions[s, 1], tce_lat.positions[s, 2]
            j = round(Int, 2 * y / tce_sq3)
            (round(Int, x - j / 2), j)
        end
        probe16.components[1] .= [iseven(i) && iseven(j) for (i, j) in tce_ij]
        @test count(probe16.components[1]) == 4
        @test order_parameter_p2x2(probe16) ≈ sqrt(3) / 4 atol = 1e-12
        for (m, n) in [(2, -1), (0, 2), (2, 1)]
            @test bragg_amplitude(probe16, m, n) == 1 / 4
        end

        # Exact grand references and the seeded igref run with the M-point
        # order parameter recorded through the hook. Tolerances: three-seed
        # calibrated (seeds 5301/5302/5303), per stat, >= 3x the maximum
        # deviation (lnXi 0.279, N 0.508, U 0.021, Psi 0.0142); sigma =
        # sqrt(16 ln 2 / 100) ~ 0.33 for lnXi. Both mu = +-0.1 corners of
        # the 800 K row sit deliberately at the |ln(z/z0)| ~ 1.45
        # trust-window edge; the low-fugacity corner (mu = -0.1,
        # z = e^{beta mu} ~ 0.23 against z0 = 1) binds, with calibrated
        # minimum Kish N_eff 23.5, so the health floor is 15, not a
        # production-grade gate
        tce_kb = 8.617333262e-5
        tce_mus = [-0.1, 0.0, 0.1]
        tce_Ts = [800.0, 1200.0, 2000.0]
        function tce_exact(mu, T)
            beta = 1 / (tce_kb * T)
            w = [exp(beta * (mu * n16[k] - e_own16[k])) for k in 1:2^16]
            Z = sum(w)
            return (log(Z), sum(w .* n16) / Z, sum(w .* e_own16) / Z, sum(w .* psi16) / Z)
        end
        Random.seed!(5301)
        tce_walkers = [LatticeWalker(deepcopy(tce_lat), energy=0.0u"eV", iter=0)
                       for _ in 1:100]
        tce_ls = LatticeGasWalkers(tce_walkers, tce_ham; assign_energy=false)
        tce_params = IdealGasReferencedGCNSParameters(mc_steps=60,
            reference_fugacity=1.0, energy_perturbation=1e-9)
        tce_df, tce_out, _ = ideal_gas_referenced_nested_sampling(tce_ls,
            tce_params, Int64(1280), MCGrandCanonicalMoves(), obs_save;
            observables=[:psi => order_parameter_p2x2])
        obs_cleanup()
        tce_liveE = [w.energy.val for w in tce_out.walkers]
        tce_liveN = [Int(sum(w.configuration.components[1])) for w in tce_out.walkers]
        tce_livePsi = Dict(:psi => [order_parameter_p2x2(w.configuration) for w in tce_out.walkers])
        tce_st = gc_thermodynamic_stats_ideal_ref(tce_df, 16, 1.0, tce_mus, tce_Ts,
            100; ω0=101 / 100, live_emax=tce_liveE, live_numbers=tce_liveN,
            observable_cols=[:psi], live_observables=tce_livePsi)
        for (i, mu) in enumerate(tce_mus), (j, T) in enumerate(tce_Ts)
            lnXi, mN, mU, mPsi = tce_exact(mu, T)
            @test tce_st.N_eff[i, j] >= 15
            @test abs(tce_st.logXi[i, j] - lnXi) < 0.9
            @test abs(tce_st.mean_N[i, j] - mN) < 1.6
            @test abs(tce_st.mean_U[i, j] - mU) < 0.07
            @test abs(tce_st.observables[:psi][i, j] - mPsi) < 0.045
        end
    end

    # ================================================================
    @testset "end-to-end: cross-Hamiltonian reweighting closure (3x3 triangular)" begin
        # Sample the reference ladder under the pair Hamiltonian A, record
        # the trio-augmented Hamiltonian B per dead point through the
        # observables hook and the occupation vectors through the dead-point
        # callback; the substituted ledger (B's energies replacing emax and
        # the live tail) closes against B's exact grand sum over the same
        # prior-volume weights, and the substituted estimator's Kish N_eff
        # never exceeds the control's
        xh_lat = obs_tri_lattice(3, 3)
        xh_hamA = GenericLatticeHamiltonian(0.0, [0.1], u"eV")
        xh_tris = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            xh_lat, [1.0, 1.0, 1.0]; expected_count=42)
        xh_hamB = ClusterLatticeHamiltonian(GenericLatticeHamiltonian(0.0, [0.1], u"eV"),
            [ClusterInteraction(0.05u"eV", xh_tris)])

        # Independent energies for both Hamiltonians over all 2^18 configs
        xh_sq3 = sqrt(3.0)
        function xh_tdist(p, q)
            best = Inf
            for ix in -1:1, iy in -1:1
                best = min(best, hypot(p[1] - q[1] + ix * 3.0, p[2] - q[2] + iy * 3 * xh_sq3))
            end
            return best
        end
        xp = [(xh_lat.positions[s, 1], xh_lat.positions[s, 2]) for s in 1:18]
        xh_pairs = [(i, j) for i in 1:18 for j in (i+1):18
                    if abs(xh_tdist(xp[i], xp[j]) - 1.0) < 1e-9]
        eA = zeros(2^18)
        eB = zeros(2^18)
        nX = zeros(Int, 2^18)
        xocc = fill(false, 18)
        for mask in 0:(2^18-1)
            for s in 1:18
                xocc[s] = ((mask >> (s - 1)) & 1) == 1
            end
            pa = 0.1 * count(pr -> xocc[pr[1]] && xocc[pr[2]], xh_pairs)
            tri = 0.05 * count(t -> xocc[t[1]] && xocc[t[2]] && xocc[t[3]], xh_tris)
            eA[mask+1] = pa
            eB[mask+1] = pa + tri
            nX[mask+1] = count_ones(mask)
        end
        xprobe = deepcopy(xh_lat)
        Random.seed!(5400)
        for _ in 1:50
            mask = rand(0:(2^18-1))
            for s in 1:18
                xprobe.components[1][s] = ((mask >> (s - 1)) & 1) == 1
            end
            @test abs(ustrip(u"eV", interacting_energy(xprobe, xh_hamA)) - eA[mask+1]) <= 1e-12
            @test abs(ustrip(u"eV", interacting_energy(xprobe, xh_hamB)) - eB[mask+1]) <= 1e-12
        end
        xh_kb = 8.617333262e-5
        function xh_exact(evec, mu, T)
            beta = 1 / (xh_kb * T)
            mx = maximum(beta * (mu * nX[k] - evec[k]) for k in 1:2^18)
            return log(sum(exp(beta * (mu * nX[k] - evec[k]) - mx) for k in 1:2^18)) + mx
        end
        xh_mus = [-0.05, 0.0, 0.05]
        xh_Ts = [800.0, 1200.0, 2000.0]

        # One seeded ladder under A. Tolerances: three-seed calibrated
        # (seeds 5401/5402/5403), per route, >= 3x the maximum deviation
        # (control 0.324, substituted 0.325); sigma = sqrt(18 ln 2 / 200)
        # ~ 0.25; the N_eff monotonicity held on every calibration seed
        Random.seed!(5401)
        xh_walkers = [LatticeWalker(deepcopy(xh_lat), energy=0.0u"eV", iter=0)
                      for _ in 1:200]
        xh_ls = LatticeGasWalkers(xh_walkers, xh_hamA; assign_energy=false)
        xh_params = IdealGasReferencedGCNSParameters(mc_steps=100,
            reference_fugacity=1.0, energy_perturbation=1e-9)
        xh_configs = Vector{Bool}[]
        xh_df, xh_out, _ = ideal_gas_referenced_nested_sampling(xh_ls, xh_params,
            Int64(2870), MCGrandCanonicalMoves(), obs_save;
            observables=[:U_b => (cfg -> ustrip(u"eV", interacting_energy(cfg, xh_hamB)))],
            dead_point_callback=(iter, w) ->
                push!(xh_configs, copy(w.configuration.components[1])))
        obs_cleanup()
        @test nrow(xh_df) > 0
        @test length(xh_configs) == nrow(xh_df)
        # The callback-collected configurations re-evaluate to the recorded
        # H_b column exactly
        xrecomp = [begin
                       xprobe.components[1] .= c
                       ustrip(u"eV", interacting_energy(xprobe, xh_hamB))
                   end for c in xh_configs]
        @test xrecomp == xh_df.U_b

        xh_liveE = [w.energy.val for w in xh_out.walkers]
        xh_liveN = [Int(sum(w.configuration.components[1])) for w in xh_out.walkers]
        xh_stA = gc_thermodynamic_stats_ideal_ref(xh_df, 18, 1.0, xh_mus, xh_Ts,
            200; ω0=201 / 200, live_emax=xh_liveE, live_numbers=xh_liveN)
        xh_dfB = copy(xh_df)
        xh_dfB.emax = copy(xh_df.U_b)
        xh_liveEB = [begin
                         xprobe.components[1] .= w.configuration.components[1]
                         ustrip(u"eV", interacting_energy(xprobe, xh_hamB))
                     end for w in xh_out.walkers]
        xh_stB = gc_thermodynamic_stats_ideal_ref(xh_dfB, 18, 1.0, xh_mus, xh_Ts,
            200; ω0=201 / 200, live_emax=xh_liveEB, live_numbers=xh_liveN)
        for (i, mu) in enumerate(xh_mus), (j, T) in enumerate(xh_Ts)
            @test abs(xh_stA.logXi[i, j] - xh_exact(eA, mu, T)) < 1.0
            @test abs(xh_stB.logXi[i, j] - xh_exact(eB, mu, T)) < 1.0
            @test 0 < xh_stB.N_eff[i, j] <= xh_stA.N_eff[i, j]
        end
    end

    # ================================================================
    @testset "end-to-end: honeycomb recipe through the generic lattice" begin
        # The honeycomb adlattice through the MLattice{1,GenericLattice}
        # inner constructor: rectangular doubled cell (orthogonal vectors
        # keep the per-axis minimum-image reduction exact) with the
        # four-site honeycomb basis. With the underlying triangular lattice
        # constant a = 1 the honeycomb edge is d_nn = a/sqrt(3) ~ 0.577 and
        # the successive neighbor distances are d_nn * (1, sqrt(3), 2,
        # sqrt(7)) = (0.577, 1.000, 1.155, 1.528) with multiplicities
        # (3, 6, 3, 6), so the ladder [0.7, 1.05, 1.30, 1.60] brackets
        # exactly one distance class per shell. Geometric cluster moves have
        # no GenericLattice method, so clusters_freq stays 0 throughout (a
        # run with clusters_freq > 0 dies mid-run with the raw MethodError;
        # documented exclusion, not exercised here).
        hc_sq3 = sqrt(3.0)
        hc_lv = [1.0 0.0 0.0; 0.0 hc_sq3 0.0; 0.0 0.0 1.0]
        hc_basis = [(0.0, 0.0, 0.0),
                    (0.5, 1 / (2 * hc_sq3), 0.0),
                    (0.5, hc_sq3 / 2, 0.0),
                    (0.0, hc_sq3 / 2 + 1 / (2 * hc_sq3), 0.0)]
        hc = MLattice{1,GenericLattice}(hc_lv, hc_basis, (6, 4, 1),
            (true, true, false), [0.7, 1.05, 1.30, 1.60],
            [fill(false, 96)], fill(true, 96))
        @test hc isa GLattice{1}
        @test num_sites(hc) == 96
        @test all(nb -> length.(nb) == [3, 6, 3, 6], hc.neighbors)

        # Bent NN-trio motif: two nearest neighbors of a common center at
        # 120 degrees, distance multiset [1/sqrt(3), 1/sqrt(3), 1]; exactly
        # 3M embeddings (three per center), checked against an independent
        # count from torus distances
        hc_dnn = 1 / hc_sq3
        function hc_tdist(p, q, Lx, Ly)
            best = Inf
            for ix in -1:1, iy in -1:1
                best = min(best, hypot(p[1] - q[1] + ix * Lx, p[2] - q[2] + iy * Ly))
            end
            return best
        end
        hp = [(hc.positions[s, 1], hc.positions[s, 2]) for s in 1:96]
        hc_own = NTuple{3,Int}[]
        for a in 1:96
            nbs = [b for b in 1:96
                   if b != a && abs(hc_tdist(hp[a], hp[b], 6.0, 4 * hc_sq3) - hc_dnn) < 1e-9]
            for x in 1:length(nbs), y in (x+1):length(nbs)
                b, c = nbs[x], nbs[y]
                abs(hc_tdist(hp[b], hp[c], 6.0, 4 * hc_sq3) - 1.0) < 1e-9 &&
                    push!(hc_own, tuple(sort([a, b, c])...))
            end
        end
        sort!(hc_own)
        @test length(hc_own) == 3 * 96
        hc_tris = enumerate_motif_embeddings(hc, [hc_dnn, hc_dnn, 1.0];
            expected_count=length(hc_own))
        @test sort(collect(hc_tris)) == hc_own

        # Full-occupancy closed form: on-site plus half the coordination-
        # weighted couplings
        hc_ham = GenericLatticeHamiltonian(-0.04,
            [-0.01, -0.005, -0.0025, -0.00125], u"eV")
        hc.components[1] .= true
        hc_closed = 96 * (-0.04) +
                    96 / 2 * sum([3, 6, 3, 6] .* [-0.01, -0.005, -0.0025, -0.00125])
        @test abs(ustrip(u"eV", interacting_energy(hc, hc_ham)) - hc_closed) < 1e-12

        # Small-cell enumeration closure: the C(16, 4) spectrum on the
        # (2, 2, 1) cell against an independent NN-bond evaluation
        hc2 = MLattice{1,GenericLattice}(hc_lv, hc_basis, (2, 2, 1),
            (true, true, false), [0.7], [vcat(fill(true, 4), fill(false, 12))],
            fill(true, 16))
        hc2_ham = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
        hc2_df, _ = exact_enumeration(hc2, hc2_ham)
        @test nrow(hc2_df) == binomial(16, 4)
        hp2 = [(hc2.positions[s, 1], hc2.positions[s, 2]) for s in 1:16]
        hc2_pairs = [(i, j) for i in 1:16 for j in (i+1):16
                     if abs(hc_tdist(hp2[i], hp2[j], 2.0, 2 * hc_sq3) - hc_dnn) < 1e-9]
        @test length(hc2_pairs) == 24
        hc2_brute = Float64[]
        hocc = fill(false, 16)
        for mask in 0:(2^16-1)
            count_ones(mask) == 4 || continue
            for s in 1:16
                hocc[s] = ((mask >> (s - 1)) & 1) == 1
            end
            push!(hc2_brute,
                  4 * (-0.04) - 0.01 * count(pr -> hocc[pr[1]] && hocc[pr[2]], hc2_pairs))
        end
        @test sort([ustrip(u"eV", e) for e in hc2_df.energy]) ≈ sort(hc2_brute) atol = 1e-12

        # A short seeded grand-canonical run on the (3, 2, 1) cell returns a
        # well-formed, monotone ledger (swap-only decorrelation)
        Random.seed!(5501)
        hc3 = MLattice{1,GenericLattice}(hc_lv, hc_basis, (3, 2, 1),
            (true, true, false), [0.7, 1.05, 1.2], [fill(false, 24)],
            fill(true, 24))
        @test all(nb -> length.(nb) == [3, 6, 3], hc3.neighbors)
        hc3_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.005, -0.0025], u"eV")
        hc3_walkers = [LatticeWalker(deepcopy(hc3), energy=0.0u"eV", iter=0)
                       for _ in 1:10]
        hc3_ls = LatticeGasWalkers(hc3_walkers, hc3_ham; assign_energy=false)
        hc3_gc = GrandCanonicalNestedSamplingParameters(mc_steps=50,
            chemical_potential=-0.05, energy_perturbation=1e-9)
        hc3_df, _, _ = grand_canonical_nested_sampling(hc3_ls, hc3_gc, Int64(300),
            MCGrandCanonicalMoves(), obs_save)
        obs_cleanup()
        @test nrow(hc3_df) > 0
        @test all(isfinite, hc3_df.omega)
        @test issorted(hc3_df.omega, rev=true)
        # Stats leg: the Omega-sorted ledger post-processes to finite
        # thermodynamics (the first GenericLattice trip through the stats)
        hc3_beta = 1 / (8.617333262e-5 * 600.0)
        hc3_E, hc3_Cv, hc3_N = gc_thermodynamic_stats(hc3_df, [hc3_beta], 10, -0.05)
        @test all(isfinite, hc3_E)
        @test all(isfinite, hc3_Cv)
        @test all(isfinite, hc3_N)
        @test 0 <= hc3_N[1] <= 24
    end
end
