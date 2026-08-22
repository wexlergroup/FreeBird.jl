@testset "SamplingSchemes Tests" begin
    
    # test exact enumeration
    include("test-exact_enumeration.jl")
    # test nested sampling
    include("test-nested_sampling.jl")
    # test NVT Monte Carlo
    include("test-nvt_monte_carlo.jl")
include("test-muvt-monte-carlo.jl")
include("test-atomistic-gc-interacting-regressions.jl")
    # test Wang-Landau
    include("test-wang_landau.jl")
    # test grand-canonical nested sampling
    include("test-grand_canonical_ns.jl")
    # test ideal-gas-referenced grand-canonical nested sampling
    include("test-ideal-gas-ref-gcns.jl")
    # seeded trajectory pins for the lattice grand-canonical kernels/drivers
    include("test-lattice-gc-trajectory-pins.jl")
    # capstone regressions for the copy-free/incremental/ledger round
    include("test-lattice-gc-scaling-regressions.jl")
    # test atomistic ideal-gas-referenced grand-canonical nested sampling
    include("test-atomistic-igref-gcns.jl")
    # absolute seeded pins for the atomistic ideal-gas-referenced driver
    include("test-atomistic-igref-driver-pins.jl")
    # test the atomistic ideal-gas-referenced ledger assembly
    include("test-atomistic-igref-stats.jl")
    # test the ideal-gas-referenced effective-sample-size reductions
    include("test-igref-ess.jl")
    # regression coverage for the continuous-space grand-canonical route
    include("test-atomistic-gc-regressions.jl")
    # test atomistic GC-NS fixed-N post-processing (ideal gas reference)
    include("test-atomistic-gcns-fixed-n.jl")
    # test lattice GC-NS fixed-N post-processing
    include("test-lattice-gcns-fixed-n.jl")
    # test hard-core lattice models (finite-J athermal recipe)
    include("test-hard-core-lattices.jl")
    # test per-dead-point observable recording and its post-processing
    include("test-ns-observables.jl")
    # test plateau-aware tie eviction in the serial atomistic NS steps
    include("test-ns-plateau-ties.jl")

end