using Test
using MembraneElectrostatics

@testset "MembraneElectrostatics.jl" begin

    
    # Test basic potential calculations
    z = 1nm
    v = V(z)
    @test typeof(v) == Float64
    @test !isnan(v)
    @test !isinf(v)
end 