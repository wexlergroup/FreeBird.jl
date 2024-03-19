@testset "Lennard-Jones Potential" begin
    @test LJParameters() == LJParameters(1.0u"eV", 1.0u"Å", Inf, 0.0u"eV")
    @test LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false) == LJParameters(0.1u"eV", 2.5u"Å", 3.5, 0.0u"eV")

    @test lj_energy(1.0u"eV", 1.0u"Å", 1.0u"Å") ≈ 0.0u"eV"
    @test lj_energy(1.0u"eV", 1.0u"Å", 2.0u"Å") ≈ -0.0615234375u"eV"
    @test lj_energy(1.0u"eV", 1.0u"Å", 0.5u"Å") ≈ 16128.0u"eV"
    @test lj_energy(1.0u"eV", 1.0u"Å", 0.0u"Å") |> isnan == true

    lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0)

    @test lj_energy(2.5u"Å", lj) ≈ -lj.shift
    @test lj_energy(0.1u"eV", 2.5u"Å", 10.0u"Å") ≈ lj.shift
    @test lj_energy(10.0u"Å", lj) ≈ 0.0u"eV" # at cutoff
    @test lj_energy(20.0u"Å", lj) ≈ 0.0u"eV" # outside cutoff

    lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=4.0, shift=0.5)

    @test lj_energy(2.5u"Å", lj) ≈ -0.5u"eV"
end