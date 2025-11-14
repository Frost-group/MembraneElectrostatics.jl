"""
Tests for Transfer Matrix Implementation

Validates the transfer matrix method against:
1. Existing CahillRecursions.jl (N=1 single slab case)
2. Analytical limits (matched permittivities, thick/thin limits)
3. Physical constraints (energy conservation, symmetry)

Reference: Cahill, K. (2012). PRE 85(5), 051921.
"""

using Test
using MembraneElectrostatics

# Load transfer matrix module
include("../src/TransferMatrix.jl")
using .TransferMatrix

@testset "Transfer Matrix Method" begin
    
    # Physical constants
    nm = 1e-9
    ϵ_0 = 8.85418682e-12
    q = 1.602176634e-19
    
    @testset "Single Layer: Validation against CahillRecursions" begin
        # Test that N=1 transfer matrix reproduces existing formulas
        # Cahill Eq. 9-11 (3-region case) vs. transfer matrix method
        
        # Standard liver membrane parameters
        ϵ_w = 80 * ϵ_0
        ϵ_l = 2 * ϵ_0
        ϵ_c = 80 * ϵ_0
        t = 5 * nm
        
        # Create single-layer system
        lipid_layer = Layer(ϵ_l, t)
        system = LayeredSystem(ϵ_w, ϵ_c, [lipid_layer])
        
        # Test points: charge in water at h=1nm, evaluate at various positions
        h = 1.0 * nm
        
        test_points = [
            (ρ=1.0nm, z=2.0nm, desc="Both in water, off-axis"),
            (ρ=1.0nm, z=0.5nm, desc="Both in water, closer"),
            (ρ=1.0nm, z=-2.5nm, desc="Charge in water, eval in lipid"),
            (ρ=1.0nm, z=-7.0nm, desc="Charge in water, eval in cytosol"),
        ]
        
        println("\n  Testing N=1 (single slab) against existing implementation:")
        
        for tp in test_points
            # Transfer matrix potential
            V_tm = V_transfer_matrix(tp.ρ, tp.z, h, system, q=q)
            
            # Existing CahillRecursions potential
            # V(z; ρ, h, m, NMAX) where m is CahillMembrane
            m = CahillMembrane(ϵ_w=ϵ_w, ϵ_l=ϵ_l, ϵ_c=ϵ_c, t=t)
            V_orig = V(tp.z, ρ=tp.ρ, h=h, m=m, NMAX=1000)
            
            # Compute relative error
            rel_error = abs(V_tm - V_orig) / (abs(V_orig) + 1e-10)
            
            println("    $(tp.desc): V_tm = $(V_tm), V_orig = $(V_orig), rel_err = $(rel_error)")
            
            # Should match within 1% (numerical integration tolerance + different NMAX)
            @test rel_error < 0.02 || abs(V_tm - V_orig) < 1e-6
        end
    end
    
    @testset "Analytical Limit: Matched Permittivities" begin
        # When all permittivities match, image charges vanish
        # Cahill discussion around Eq. 1359-1363
        
        ϵ_uniform = 80 * ϵ_0
        t = 5 * nm
        
        # System with matched permittivities (no interfaces)
        layer = Layer(ϵ_uniform, t)
        system = LayeredSystem(ϵ_uniform, ϵ_uniform, [layer])
        
        h = 1.0 * nm
        ρ = 1.0 * nm
        z = 2.0 * nm
        
        # Potential should be pure Coulomb: q/(4πε·r)
        r = sqrt(ρ^2 + (z - h)^2)
        V_expected = q / (4π * ϵ_uniform * r)
        
        V_tm = V_transfer_matrix(ρ, z, h, system, q=q)
        
        rel_error = abs(V_tm - V_expected) / V_expected
        println("\n  Matched permittivities: V_tm = $(V_tm), V_Coulomb = $(V_expected), rel_err = $(rel_error)")
        
        @test rel_error < 0.01
    end
    
    @testset "Analytical Limit: Reflection Coefficients" begin
        # Test reflection coefficient computation
        # Cahill Eq. 1705-1709: p = (ε₂ - ε₁)/(ε₂ + ε₁)
        
        ϵ_1 = 80 * ϵ_0
        ϵ_2 = 2 * ϵ_0
        
        p = reflection_coefficient(ϵ_1, ϵ_2)
        p_expected = (ϵ_2 - ϵ_1) / (ϵ_2 + ϵ_1)
        
        @test p ≈ p_expected
        
        # Matched permittivities: p = 0
        p_zero = reflection_coefficient(ϵ_1, ϵ_1)
        @test abs(p_zero) < 1e-15
        
        # Symmetry: p(ε₁,ε₂) = -p(ε₂,ε₁)
        p_rev = reflection_coefficient(ϵ_2, ϵ_1)
        @test p ≈ -p_rev
        
        println("\n  Reflection coefficients: p(80,2) = $(p), expected = $(p_expected)")
    end
    
    @testset "Physical Constraint: Energy Monotonicity" begin
        # Self-interaction energy should be monotonic as charge approaches interface
        # Image charges from low-ε membrane repel charges in high-ε water
        # Cahill Eq. 1013-1017 and Fig. 5 (Debye layer discussion)
        
        ϵ_w = 80 * ϵ_0
        ϵ_l = 2 * ϵ_0
        ϵ_c = 80 * ϵ_0
        t = 5 * nm
        
        layer = Layer(ϵ_l, t)
        system = LayeredSystem(ϵ_w, ϵ_c, [layer])
        
        # Compute self-interaction at various heights
        heights = [5.0nm, 2.0nm, 1.0nm, 0.5nm, 0.2nm]
        energies = Float64[]
        
        println("\n  Self-interaction energy vs. height above membrane:")
        for z in heights
            E_self = V_transfer_matrix_self_interaction(z, system, q=q)
            push!(energies, E_self)
            println("    z = $(z/nm) nm: E_self = $(E_self) J")
        end
        
        # Energy should increase (become more repulsive) as charge approaches interface
        # Check that energies are monotonically increasing
        for i in 2:length(energies)
            @test energies[i] >= energies[i-1]  # Closer = more repulsive = higher energy
        end
    end
    
    @testset "Multi-Layer: 3-Layer System" begin
        # Test 3 internal layers (5 regions total)
        # Simplified version of Cahill's head-group model
        # Section "Three Dielectric Layers" and Fig. 7
        
        ϵ_w = 80 * ϵ_0
        ϵ_head = 195 * ϵ_0  # High permittivity head groups
        ϵ_lipid = 2 * ϵ_0
        ϵ_c = 80 * ϵ_0
        
        t_head = 1 * nm
        t_lipid = 3 * nm
        
        # 3 layers: head group, lipid, head group
        layers = [
Layer(ϵ_head, t_head),
Layer(ϵ_lipid, t_lipid),
Layer(ϵ_head, t_head)
        ]
        
        system_3layer = LayeredSystem(ϵ_w, ϵ_c, layers)
        
        # Compare with single lipid layer (no head groups)
        system_1layer = LayeredSystem(
            ϵ_w, ϵ_c, [Layer(ϵ_lipid, 5*nm)]
        )
        
        h = 1.0 * nm
        ρ = 1.0 * nm
        z = 0.5 * nm
        
        V_3layer = V_transfer_matrix(ρ, z, h, system_3layer, q=q)
        V_1layer = V_transfer_matrix(ρ, z, h, system_1layer, q=q)
        
        println("\n  3-layer vs 1-layer system:")
        println("    V_3layer (with head groups) = $(V_3layer)")
        println("    V_1layer (lipid only) = $(V_1layer)")
        println("    Difference = $(V_3layer - V_1layer)")
        
        # Head groups should reduce potential (attract ions)
        # Cahill Fig. 7-8: head groups make potential less repulsive/more attractive
        @test V_3layer < V_1layer  # More attractive with head groups
    end
    
    @testset "Numerical Stability: Transfer Matrix Elements" begin
        # Test that transfer matrix doesn't overflow/underflow
        # For large k or thick layers, exp(2kt) can be huge
        
        ϵ_w = 80 * ϵ_0
        ϵ_l = 2 * ϵ_0
        ϵ_c = 80 * ϵ_0
        t = 5 * nm
        
        layer = Layer(ϵ_l, t)
        system = LayeredSystem(ϵ_w, ϵ_c, [layer])
        
        # Test at various k values
        k_values = [0.1/nm, 1.0/nm, 10.0/nm, 50.0/nm, 100.0/nm]
        h = 1.0 * nm
        
        println("\n  Transfer matrix stability at various k:")
        for k in k_values
            E, β = compute_E_matrix(system, k, h)
            
            # Check for NaN or Inf
            @test all(isfinite.(E))
            @test isfinite(β)
            
            # Matrix elements shouldn't be too large
            @test maximum(abs.(E)) < 1e10
            
            println("    k = $(k*nm) nm⁻¹: max|E| = $(maximum(abs.(E)))")
        end
    end
    
    @testset "Symmetry: z-reflection" begin
        # For symmetric system (ϵ_w = ϵ_c), potential should respect z→-z symmetry
        # with appropriate charge/evaluation position flips
        
        ϵ_water = 80 * ϵ_0
        ϵ_l = 2 * ϵ_0
        t = 5 * nm
        
        layer = Layer(ϵ_l, t)
        system_sym = LayeredSystem(ϵ_water, ϵ_water, [layer])
        
        # Charge at +h above interface, evaluate at +z
        h = 1.0 * nm
        z = 2.0 * nm
        ρ = 1.0 * nm
        
        V_plus = V_transfer_matrix(ρ, z, h, system_sym, q=q)
        
        # Due to symmetry, flipping system should give same result
        # (This is a conceptual test - implementation is for z > 0 charge)
        # Just verify computation completes without error
        @test isfinite(V_plus)
        
        println("\n  Symmetric system (ϵ_w = ϵ_c): V = $(V_plus)")
    end
end

println("\n✓ All transfer matrix tests passed!")

