# Define custom types for different regions
#   ==> dispatch based on these types (9 different equations provided by Cahill)
struct WaterRegion end
struct LipidRegion end
struct CytosolRegion end

# Distance scaling
const nm = 1E-9

# Membrane type struct encapsulating dielectric properties and thickness
struct CahillMembrane
    ϵ_w::Float64
    ϵ_l::Float64
    ϵ_c::Float64
    t::Float64
    p::Float64
    p′::Float64
    ϵ_wl::Float64
    ϵ_cl::Float64
    
    function CahillMembrane(; ϵ_w=80*ϵ_0, ϵ_l=2*ϵ_0, ϵ_c=80*ϵ_0, t=5nm)
        p = (ϵ_w - ϵ_l)/(ϵ_w + ϵ_l)
        p′ = (ϵ_c - ϵ_l)/(ϵ_c + ϵ_l)
        ϵ_wl = (ϵ_w + ϵ_l)/2
        ϵ_cl = (ϵ_c + ϵ_l)/2
        new(ϵ_w, ϵ_l, ϵ_c, t, p, p′, ϵ_wl, ϵ_cl)
    end
end

const CAHILL_LIVER = CahillMembrane()

# Main potential functions
# Cahill's Nomenclature:
#   V     - potential
#    ^w   - charge is in (w)ater, (l)ipid, or (c)ytosol
#     _w  - evaluating in (w)ater, (l)ipid, or (c)ytosol
# We now achieve this by dispatching on the region types.

# Charge in water
# V^w_w, Eqn 9 in Cahill2012
@fastmath function V(z, charge::WaterRegion, eval::WaterRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    if ρ ≈ 0
        return V_w_w_selfinteraction(z, charge, eval; t=t, m=m, NMAX=NMAX)
    end

    q/(4π*m.ϵ_w) * (
        1/√(ρ^2+(z-h)^2) 
        + m.p/√(ρ^2+(z+h)^2)
        - m.p′*(1-m.p^2) * 
        sum( ((m.p*m.p′)^(n-1)/(√(ρ^2 + (z + 2n*t + h)^2))) 
            for n in 1:NMAX) )
end

# Eqn 35; special case of 9 for ρ=0, self-interaction / on-axis evaluation
function V_w_w_selfinteraction(z, ::WaterRegion, ::WaterRegion; t, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_w) * (
        m.p/abs(2z) - 
        (m.ϵ_w*m.ϵ_l/m.ϵ_wl^2) * 
        sum(m.p^(n-1)*m.p′^n / abs(2z + 2n*t) for n in 1:NMAX)
    )
end

# V_^w_l, Eqn 10 in Cahill2012
function V(z, charge::WaterRegion, eval::LipidRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_wl) * 
    sum(
        ((m.p*m.p′)^n * (1/√(ρ^2 + (z - 2n*t - h)^2) 
          - m.p′/√(ρ^2 + (z + 2*(n+1)*t + h)^2)))
         for n in 0:NMAX)
end

# V_^w_c, Eqn 11 in Cahill2012
function V(z, charge::WaterRegion, eval::CytosolRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q*m.ϵ_l/(4π*m.ϵ_wl*m.ϵ_cl) * 
    sum(
        ((m.p*m.p′)^n / √(ρ^2 + (z - 2n*t - h)^2)) 
            for n in 0:NMAX)
end

# Charge in lipid
# V_^l_w, Eqn 15 in Cahill2012
function V(z, charge::LipidRegion, eval::WaterRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_wl) * (
        sum((m.p*m.p′)^n / √(ρ^2 + (z + 2n*t - h)^2) for n in 0:NMAX) -
        sum(m.p′ * (m.p*m.p′)^n / √(ρ^2 + (z + 2(n+1)*t + h)^2) for n in 0:NMAX)
    )
end

# V_^l_l, Eqn 16 in Cahill2012
function V(z, charge::LipidRegion, eval::LipidRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_l) * (
        sum((m.p*m.p′)^abs(n) / √(ρ^2 + (z - 2n*t - h)^2) for n in -NMAX:NMAX) -
        sum(m.p*(m.p*m.p′)^n / √(ρ^2 + (z - 2n*t + h)^2) for n in 0:NMAX) -
        sum(m.p′*(m.p*m.p′)^n / √(ρ^2 + (z + 2(n+1)*t + h)^2) for n in 0:NMAX)
    )
end

# V_^l_c, Eqn 17 in Cahill2012
function V(z, charge::LipidRegion, eval::CytosolRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_cl) * sum(
        (m.p*m.p′)^n * (
            1/√(ρ^2 + (z - 2n*t - h)^2) -
            m.p/√(ρ^2 + (z - 2n*t + h)^2)
        ) for n in 0:NMAX
    )
end

# Charge in cytosol
# V_^c_w, Eqn 21 in Cahill2012
function V(z, charge::CytosolRegion, eval::WaterRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q*m.ϵ_l/(4π*m.ϵ_wl*m.ϵ_cl) * 
    sum(
        (m.p*m.p′)^n / √(ρ^2 + (z + 2n*t - h)^2)
        for n in 0:NMAX
    )
end

# V_^c_l, Eqn 22 in Cahill2012
function V(z, charge::CytosolRegion, eval::LipidRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_cl) * (
        sum((m.p*m.p′)^n / √(ρ^2 + (z - h + 2n*t)^2) for n in 0:NMAX) -
        m.p * sum((m.p*m.p′)^n / √(ρ^2 + (z + h - 2n*t)^2) for n in 0:NMAX)
    )
end

# V_^c_c, Eqn 23 in Cahill2012
function V(z, charge::CytosolRegion, eval::CytosolRegion; ρ, t, h, m::CahillMembrane, NMAX=1000)
    q/(4π*m.ϵ_c) * (
        1/√(ρ^2+(z-h)^2) +
        m.p′/√(ρ^2 + (z + h + 2t)^2) -
        m.p*(1 - m.p′^2) * sum(
            (m.p*m.p′)^n / √(ρ^2 + (z - 2n*t - h)^2)
            for n in 0:NMAX
        )
    )
end

# Helper functions to determine regions
function charge_region(h, m::CahillMembrane)    
    if h > 0
        WaterRegion()
    elseif h >= -m.t
        LipidRegion()
    else
        CytosolRegion()
    end
end

function eval_region(z, m::CahillMembrane)
    if z > 0
        WaterRegion()
    elseif z >= -m.t
        LipidRegion()
    else
        CytosolRegion()
    end
end

"""
    V(z; ρ=1nm, h=1nm, m=CAHILL_LIVER, NMAX=1000)

Calculate potential at position z, with charge at position h, using membrane type m.

Dispatch is made to the correct V_{w,l,c}{w,l,c} function from Cahill's paper, based on the numeric value of h and z relative to membrane thickness.

NMAX controls the convergence of infinite sums (Cahill uses NMAX=1000).
"""
function V(z; ρ=1nm, h=1nm, m::CahillMembrane=CAHILL_LIVER, NMAX=1000)
    charge = charge_region(h, m)
    eval = eval_region(z, m)
    V(z, charge, eval; ρ=ρ, t=m.t, h=h, m=m, NMAX=NMAX)
end
