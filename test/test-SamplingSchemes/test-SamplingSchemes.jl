@testset "SamplingSchemes Tests" begin
    
    # test exact enumeration
    include("test-exact_enumeration.jl")
    # test nested sampling
    include("test-nested_sampling.jl")
    # test NVT Monte Carlo
    include("test-nvt_monte_carlo.jl")
    # test Wang-Landau
    include("test-wang_landau.jl")
    # test grand-canonical nested sampling
    include("test-grand_canonical_ns.jl")
    # test ideal-gas-referenced grand-canonical nested sampling
    include("test-ideal-gas-ref-gcns.jl")

end