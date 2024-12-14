@testset "Lennard-Jones Potential" begin
    @test LJParameters() == LJParameters(1.0u"eV", 1.0u"Å", Inf, 0.0u"eV")
    @test typeof(LJParameters().epsilon) == typeof(1.0u"eV")
    @test typeof(LJParameters().sigma) == typeof(1.0u"Å")
    @test typeof(LJParameters().cutoff) == Float64
    @test typeof(LJParameters().shift) == typeof(0.0u"eV")

    @test typeof(LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=4.5)) == LJParameters
    @test LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false) == LJParameters(0.1u"eV", 2.5u"Å", 3.5, 0.0u"eV")

    @test lj_energy(1.0u"eV", 1.0u"Å", 1.0u"Å") == 0.0u"eV"
    @test lj_energy(1.0u"eV", 1.0u"Å", 2.0u"Å") ≈ -0.0615234375u"eV"
    @test lj_energy(1.0u"eV", 1.0u"Å", 0.5u"Å") ≈ 16128.0u"eV"
    @test lj_energy(1.0u"eV", 1.0u"Å", 0.0u"Å") |> isnan == true

    lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0)
    @test lj.shift == lj_energy(0.1u"eV", 2.5u"Å", 10.0u"Å")

    @test lj_energy(2.5u"Å", lj) == lj_energy(0.1u"eV", 2.5u"Å", 2.5u"Å") - lj.shift   # at the minima
    @test lj_energy(2.5u"Å", lj) == -lj.shift   # at the minima
    @test lj_energy(6.0u"Å", lj) ≈ -0.00198452714u"eV"    # at cutoff point
    @test lj_energy(20.0u"Å", lj) ≈ 0.0u"eV"    # outside of cutoff

    lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0, shift=0.5)

    @test lj_energy(2.5u"Å", lj) == -0.5u"eV"
    

    v1 = [lj, lj, lj, lj, lj, lj, lj, lj, lj]
    v2 = [lj, lj, lj, lj, lj, lj]

    @test CompositeLJParameters(3, v1) |> typeof == CompositeLJParameters{3}
    @test CompositeLJParameters(3, v2) |> typeof == CompositeLJParameters{3}
    @test CompositeLJParameters(3, v1) == CompositeLJParameters(3, v2)

    
    ljs = [LJParameters(epsilon=e) for e in [0.11, 0.21, 0.31, 0.12, 0.22, 0.32, 0.13, 0.23, 0.33]]
    @test CompositeLJParameters(3, ljs) == CompositeLJParameters{3}(reshape(ljs, 3, 3))

    ljs1 = [LJParameters(epsilon=e) for e in [0.11, 0.12, 0.13, 0.22, 0.23, 0.33]]
    ljs2 = [LJParameters(epsilon=e) for e in [0.11, 0.12, 0.13, 0.12, 0.22, 0.23, 0.13, 0.23, 0.33]]

    @test CompositeLJParameters(3, ljs1) == CompositeLJParameters(3, ljs2)

end