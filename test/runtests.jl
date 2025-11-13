using Test
using MembraneElectrostatics

@testset "MembraneElectrostatics.jl" begin

    
    # Test basic potential calculations
    z = 1nm
    v = V(z)
    @test typeof(v) == Float64
    @test !isnan(v)
    @test !isinf(v)

    @testset "Figure 1 calculations" begin
        # Use same parameters as examples /figure1.jl
        m = CahillMembrane(t=5nm)
        Zs = collect(-(m.t+5nm):m.t/5:5nm) # smaller number of data points

        # Test that all potentials evaluate without NaN or Inf
        for h in [-6nm, 1nm]
            for z in Zs
                v = V(z, h=h, m=m, NMAX=10)
                @test !isnan(v)
                @test !isinf(v)
            end
        end
    end
end

@testset "IonMonteCarlo.jl" begin
    include("IonMonteCarlo.jl")
end
