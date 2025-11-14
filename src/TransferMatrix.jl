"""
Transfer Matrix Method for N-Layer Dielectric Systems

Implements the transfer matrix formalism from:
Cahill, K. (2012). Models of membrane electrostatics. Physical Review E, 85(5), 051921.
Section VIII "Several Dielectric Layers" 

Arxiv https://arxiv.org/abs/1101.4265 'Tex Source' - qbiopre.tex lines 1584-1774

This module provides automatic computation of electrostatic potentials for arbitrary
N-layer dielectric systems using the transfer matrix method described by Cahill.

The method systematically generates image charge series through matrix multiplication. 
"""
module TransferMatrix

using SpecialFunctions: besselj0
using QuadGK: quadgk

export Layer, LayeredSystem
export reflection_coefficient, transfer_matrix, compute_E_matrix
export reflection_coefficients, V_transfer_matrix, V_transfer_matrix_self_interaction
export compare_with_single_slab

"""
    Layer

Represents a dielectric layer with permittivity ϵ and thickness t.

# Fields
- `ϵ::Float64`: Electric permittivity (F/m)
- `t::Float64`: Layer thickness (m, typically nm scale)
"""
struct Layer
    ϵ::Float64
    t::Float64
end

"""
    LayeredSystem

Describes a system of N dielectric layers sandwiched between two semi-infinite regions.

The geometry is:
- z > 0: Semi-infinite region with permittivity ϵ_w ('water')
- -t₁ < z < 0: Layer 1 with permittivity ϵ₁, thickness t₁
- -t₁-t₂ < z < -t₁: Layer 2 with permittivity ϵ₂, thickness t₂
- ... (N layers total)
- z < -Σtᵢ: Semi-infinite region with permittivity ϵ_c (cytosol)

Reference: Cahill Eqn. 45 onwards

# Fields
- `ϵ_w::Float64`: Permittivity of water (semi-infinite, z > 0)
- `ϵ_c::Float64`: Permittivity of cytosol (semi-infinite, z < -Σtᵢ)
- `layers::Vector{Layer}`: Internal dielectric layers (ordered from z=0 downward)
"""
struct LayeredSystem
    ϵ_w::Float64
    ϵ_c::Float64
    layers::Vector{Layer}
end

"""
    reflection_coefficient(ϵ₁, ϵ₂)

Compute reflection coefficient at interface between two dielectrics.

Cahill Eqn. 49: pᵢ = (εᵢ₊₁ - εᵢ)/(εᵢ₊₁ + εᵢ)

It's a bit unclear but these define reflection/ transmission coefficients so you can stack them up. 

This coefficient encodes the strength of image charges at the interface.
- p > 0: Higher permittivity on far side (attractive image)
- p < 0: Lower permittivity on far side (repulsive image)
- p = 0: Matched permittivities (no reflection)

# Arguments
- `ϵ₁::Float64`: Permittivity of first medium
- `ϵ₂::Float64`: Permittivity of second medium

# Returns
- `Float64`: Reflection coefficient p ∈ (-1, 1)
"""
reflection_coefficient(ϵ₁, ϵ₂) = (ϵ₂ - ϵ₁) / (ϵ₂ + ϵ₁)

"""
    transfer_matrix(pᵢ, yᵢ, ε̄ᵢ, εᵢ₊₁)

Compute 2×2 transfer matrix for propagation across interface i.

Cahill Eqn. 50 (lines 1725-1733 in the pre-print Tex):
```
[mᵢ₊₁]     [ε̄ᵢ/εᵢ₊₁] * [1       pᵢyᵢ   ] * [mᵢ]
[fᵢ₊₁]  =              [pᵢ/yᵢ     1     ]   [fᵢ]
```

where:
- mᵢ, fᵢ are coefficients of e^(kz) and e^(-kz) in layer i
- pᵢ = (εᵢ₊₁ - εᵢ)/(εᵢ₊₁ + εᵢ) is reflection coefficient
- yᵢ = exp[2k(t₁ + ... + tᵢ)] is phase accumulation
- ε̄ᵢ = (εᵢ₊₁ + εᵢ)/2 is average permittivity

# Arguments
- `pᵢ::Float64`: Reflection coefficient at interface i
- `yᵢ::Float64`: Phase factor exp(2k × cumulative_thickness)
- `ε̄ᵢ::Float64`: Average permittivity (εᵢ₊₁ + εᵢ)/2
- `εᵢ₊₁::Float64`: Permittivity of next layer

# Returns
- `Matrix{Float64}`: 2×2 transfer matrix
"""
function transfer_matrix(pᵢ::Float64, yᵢ::Float64, ε̄ᵢ::Float64, εᵢ₊₁::Float64)
    prefactor = ε̄ᵢ / εᵢ₊₁
    return prefactor * [1.0    pᵢ*yᵢ
                        pᵢ/yᵢ  1.0]
end

"""
    compute_E_matrix(system::LayeredSystem, k::Float64, h::Float64)

Compute total transfer matrix E^(n) by multiplying matrices for all interfaces.

Cahill Eqn. 51 (lines 1747-1754 in the pre-print Tex):
```
E^(n) = ∏ᵢ₌₀ⁿ (ε̄ᵢ/εᵢ₊₁) * [1  pᵢyᵢ; pᵢ/yᵢ  1]
```

This matrix relates the boundary coefficients in 'water' (top of stack) to those in 'cytosol' (bottom of stack):
```
[mₙ₊₁]       [β]
[fₙ₊₁]  = E^(n) * [u]
```

where β = q·exp(-kh)/(4πεw) represents the source charge.

# Arguments
- `system::LayeredSystem`: The layered dielectric system
- `k::Float64`: Wavevector in Hankel transform (units: 1/m)
- `h::Float64`: Height of source charge above interface (m)

# Returns
- `Matrix{Float64}`: 2×2 product matrix E^(n)
- `Float64`: Source term β(k)
"""
function compute_E_matrix(system::LayeredSystem, k::Float64, h::Float64)
    N = length(system.layers)
    
    # Source term β = q·exp(-kh)/(4πεw)
    # We'll factor out q/(4πεw) later, so just track exp(-kh)
    β_factor = exp(-k * h)
    
    # Cumulative thickness for computing yᵢ = exp(2k·Σtⱼ)
    cumulative_t = 0.0
    
    # Build permittivity sequence: ε₀=εw, ε₁, ε₂, ..., εₙ, εₙ₊₁=εc
    ε_seq = vcat(system.ϵ_w, [layer.ϵ for layer in system.layers], system.ϵ_c)
    
    # Start with identity matrix
    E = [1.0 0.0; 0.0 1.0]
    # Multiply transfer matrices for each interface (i = 0 to n), Eqn. 51
    for i in 0:N
        εᵢ = ε_seq[i+1]      # Julia 1-based indexing
        εᵢ₊₁ = ε_seq[i+2]
        
        # Cahill Eqn. 49
        pᵢ = reflection_coefficient(εᵢ, εᵢ₊₁)
        ε̄ᵢ = (εᵢ + εᵢ₊₁) / 2
        
        # Update cumulative thickness
        if i > 0
            cumulative_t += system.layers[i].t
        end
        
        # Cahill inline, just below 50: yᵢ = exp[2k(t₁ + ... + tᵢ)]
        #  It's quite mad it doesn't get its own eqn number, as this phase factor is absolutely key
        # For i=0, y₀ = 1
        yᵢ = (i == 0) ? 1.0 : exp(2 * k * cumulative_t)
        
        # Multiply current matrix by transfer matrix for this interface
        # Cahill Eqn. 51
        Tᵢ = transfer_matrix(pᵢ, yᵢ, ε̄ᵢ, εᵢ₊₁)
        E = Tᵢ * E
    end
    
    return E, β_factor
end

"""
    reflection_coefficients(system::LayeredSystem, k::Float64, h::Float64)

Extract reflection and transmission coefficients from transfer matrix.

Cahill Eqn. 53:
```
u(k) = -β · E₂₁^(n) / E₂₂^(n)
d(k) = β · (E₁₁^(n) - E₁₂^(n)·E₂₁^(n)/E₂₂^(n))
```

where:
- u(k): Reflection coefficient (upward-propagating wave in 'water')
- d(k): Transmission coefficient (downward-propagating wave in 'cytosol')
- β(k) = q·exp(-kh)/(4πεw): Source term

These coefficients encode all image charge contributions through the layered structure.

# Arguments
- `system::LayeredSystem`: The layered dielectric system
- `k::Float64`: Wavevector (1/m)
- `h::Float64`: Source height (m)

# Returns
- `Tuple{Float64, Float64}`: (u(k), d(k)) reflection and transmission coefficients
"""
function reflection_coefficients(system::LayeredSystem, k::Float64, h::Float64)
    E, β_factor = compute_E_matrix(system, k, h)
    
    # Extract matrix elements (Cahill notation)
    E₁₁ = E[1, 1]
    E₁₂ = E[1, 2]
    E₂₁ = E[2, 1]
    E₂₂ = E[2, 2]
    
    # Cahill Eqn. 53
    # Note: β_factor already includes exp(-kh), full β = q·exp(-kh)/(4πεw)
    u = -β_factor * E₂₁ / E₂₂
    d = β_factor * (E₁₁ - E₁₂ * E₂₁ / E₂₂)
    
    return u, d
end

"""
    V_transfer_matrix(ρ, z, h, system::LayeredSystem; 
                      q=q, k_max=100/nm, n_points=500, region=:water)

Compute electrostatic potential using transfer matrix method and inverse Hankel transform.

Cahill Eqn. 45:
```
V^w_w(ρ,z) = ∫₀^∞ dk J₀(kρ) [q/(4πεw)·e^(-k|z-h|) + u(k)·e^(-kz)]
```

The potential consists of:
1. Direct Coulomb term: q/(4πεw)·1/√(ρ² + (z-h)²)
2. Image charge contributions: encoded in u(k) via transfer matrix

# Arguments
- `ρ::Float64`: Radial distance from charge (m)
- `z::Float64`: Vertical position where potential is evaluated (m)
- `h::Float64`: Vertical position of source charge (m, h > 0 for water)
- `system::LayeredSystem`: The layered dielectric system

# Keyword Arguments
- `q::Float64`: Charge magnitude (C), default from MembraneElectrostatics
- `k_max::Float64`: Maximum wavevector for integration (1/m), default 100/nm
- `n_points::Int`: Number of quadrature points, default 500
- `region::Symbol`: Region where potential is evaluated (:water, :cytosol, :layer)

# Returns
- `Float64`: Electrostatic potential (V) at position (ρ, z)

# Notes
- Integration uses adaptive quadrature for oscillatory Bessel function
- Direct term diverges as ρ,z→h; handled separately
- Numerical precision depends on k_max and n_points
"""
function V_transfer_matrix(ρ::Float64, z::Float64, h::Float64, 
                          system::LayeredSystem;
                          q::Float64=1.602176634e-19,  # electron charge
                          k_max::Float64=100.0/1e-9,   # 100/nm in 1/m
                          n_points::Int=500,
                          region::Symbol=:water)
    
    ϵ_w = system.ϵ_w
    
    # Cahill Eqn. 45: Bessel expansion of potential
    # V = ∫₀^∞ dk J₀(kρ) [β·e^(-k|z-h|) + u(k)·e^(-kz)]
    # where β = q/(4πεw) for the direct term (handled separately)
    
    # Direct Coulomb term: q/(4πεw·r) where r = √(ρ² + (z-h)²)
    # Cahill Eqn. 3: 1/r = ∫₀^∞ dk J₀(kρ) e^(-k|z-h|)
    r = sqrt(ρ^2 + (z - h)^2)
    V_direct = (r > 1e-15) ? q / (4π * ϵ_w * r) : 0.0
    
    # Image charge contribution via transfer matrix
    # Integrate: ∫₀^∞ dk J₀(kρ) u(k) e^(-kz)
    function integrand(k::Float64)
        if k < 1e-10
            return 0.0  # Avoid numerical issues at k=0
        end
        
        # For very high k, exp(-kz) provides natural cutoff
        # If exp(-kz) is negligible, skip expensive transfer matrix calculation
        exponential_term = exp(-k * abs(z))
        if exponential_term < 1e-100
            return 0.0
        end
        
        # Get reflection coefficient from transfer matrix
        u, _ = reflection_coefficients(system, k, h)
        
        # Check for numerical overflow in transfer matrix
        if !isfinite(u)
            return 0.0  # At very high k, contribution is negligible anyway
        end
        
        # Cahill Eqn. 45 (lines 1610-1613 in the pre-print Tex): J₀(kρ) u(k) e^(-kz)
        # Note: u already includes q/(4πεw)·exp(-kh) factor from β
        bessel_term = besselj0(k * ρ)
        
        result = bessel_term * u * exponential_term * q / (4π * ϵ_w)
        
        # Final safety check
        return isfinite(result) ? result : 0.0
    end
    
    # Adaptive quadrature integration
    # Bessel function oscillates; may need higher k_max for small ρ
    V_image, err = quadgk(integrand, 0.0, k_max, rtol=1e-6, atol=1e-10)
    
    return V_direct + V_image
end

"""
    V_transfer_matrix_self_interaction(z, system::LayeredSystem; 
                                       q=q, k_max=100/nm, n_points=500)

Compute self-interaction potential (image charges only, no direct 1/r term).

Cahill Eqn. 61 (lines 12011-12041 in the pre-print Tex): Self-potential felt by
charge due to induced polarization.

For a charge at height z, this computes the potential at z due to all image
charges, excluding the infinite direct term. Used in Monte Carlo energy
calculations.

# Arguments
- `z::Float64`: Position of charge (m)
- `system::LayeredSystem`: The layered dielectric system

# Keyword Arguments
- `q::Float64`: Charge magnitude (C)
- `k_max::Float64`: Maximum wavevector (1/m)
- `n_points::Int`: Quadrature points

# Returns
- `Float64`: Self-interaction energy q·V_self (J)
"""
function V_transfer_matrix_self_interaction(z::Float64, system::LayeredSystem;
                                           q::Float64=1.602176634e-19,
                                           k_max::Float64=100.0/1e-9,
                                           n_points::Int=500)
    
    # Self-interaction: charge at z, evaluate potential at z, exclude (infinite)direct term
    # This is purely the image charge contribution
    
    ϵ_w = system.ϵ_w
    
    function integrand(k::Float64)
        if k < 1e-10
            return 0.0
        end
        
        # Early cutoff for negligible contributions
        exponential_term = exp(-k * abs(z))
        if exponential_term < 1e-100
            return 0.0
        end
        
        # For self-interaction, h = z (source and evaluation at same point)
        u, _ = reflection_coefficients(system, k, z)
        
        # Handle numerical overflow
        if !isfinite(u)
            return 0.0
        end
        
        # Only the image term, no direct Coulomb
        result = u * exponential_term / (4π * ϵ_w)
        
        return isfinite(result) ? result : 0.0
    end
    
    V_image, err = quadgk(integrand, 0.0, k_max, rtol=1e-6, atol=1e-10)
    
    # Return energy: q × V
    return q * V_image
end


end # module TransferMatrix

