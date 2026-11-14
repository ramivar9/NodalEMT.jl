# Test suite of NodalEMT.jl.  Run from the repository root with
#     julia --project=. --threads=4 -e "using Pkg; Pkg.test()"
# The test systems are read from data/psse/ with PowerSystems.jl and PowerFlows.jl.

using Test

include("test_utils.jl")

@testset verbose = true "NodalEMT.jl" begin
    @testset "PSS/E cases" begin include("test_psse.jl") end
    @testset "Bergeron lines" begin include("test_bergeron.jl") end
    @testset "compensation" begin include("test_compensation.jl") end
    @testset "scenario batches" begin include("test_batch.jl") end
    @testset "model sets" begin include("test_models.jl") end
    @testset "simulation" begin include("test_simulation.jl") end
end
