# Define custom types for different regions
#   ==> dispatch based on these types (9 different equations provided by Cahill)
struct WaterRegion end
struct LipidRegion end
struct CytosolRegion end

# Material constants
"Water absolute dielectric."
const ϵ_w = 80 * ϵ_0
"Lipid absolute dielectric."
const ϵ_l = 2 * ϵ_0
"Cytosol absolute dielectric."
const ϵ_c = 80 * ϵ_0

# Distance scaling
const nm = 1E-9

# Helper constants
const p = (ϵ_w - ϵ_l)/(ϵ_w + ϵ_l)
const p′ = (ϵ_c - ϵ_l)/(ϵ_c + ϵ_l)
const ϵ_wl = (ϵ_w + ϵ_l)/2
const ϵ_cl = (ϵ_c + ϵ_l)/2

# Main potential functions
# Cahill's Nomenclature:
#   V_   - potential
#     w  - charge is in (w)ater, (l)ipid, or (c)ytosol
#     w  - evaluating in (w)ater, (l)ipid, or (c)ytosol
# We now achieve this by dispatching on the region types.
"""
    V(z, charge::T, eval::U; ρ, t, h, NMAX) where {T,U}

Calculate potential for charge in region T, evaluating in region U.
NMAX controls the convergence of infinite sums.
"""
function V(z, charge::WaterRegion, eval::WaterRegion; ρ, t, h, NMAX=1000)
    q/(4π*ϵ_w) * (
        1/√(ρ^2+(z-h)^2) 
        + p/√(ρ^2+(z+h)^2)
        - p′*(1-p^2) * 
        sum( (p*p′^(n-1)/(√(ρ^2 + (z + 2n*t + h)^2))) 
            for n in 1:NMAX) )
end

function V(z, charge::WaterRegion, eval::LipidRegion; ρ, t, h, NMAX=1000)
    q/(4π*ϵ_wl) * 
    sum(
        ((p*p′)^n * (1/√(ρ^2 + (z - 2n*t - h)^2) 
          - p′/√(ρ^2 + (z + 2*(n+1)*t + h)^2)))
         for n in 0:NMAX)
end

function V(z, charge::WaterRegion, eval::CytosolRegion; ρ, t, h, NMAX=1000)
    q*ϵ_l/(4π*ϵ_wl*ϵ_cl) * 
    sum(
        ((p*p′)^n / √(ρ^2 + (z - 2n*t - h)^2)) 
            for n in 0:NMAX)
end

# Add methods for charge in Lipid and Cytosol...
# function V(z, charge::LipidRegion, eval::WaterRegion; kwargs...)
# function V(z, charge::LipidRegion, eval::LipidRegion; kwargs...)
# function V(z, charge::LipidRegion, eval::CytosolRegion; kwargs...)
# function V(z, charge::CytosolRegion, eval::WaterRegion; kwargs...)
# etc...

# Helper functions to determine regions
function charge_region(h, t)    
    if h > 0
        WaterRegion()
    elseif h > -t
        LipidRegion()
    else
        CytosolRegion()
    end
end

function eval_region(z, t)
    if z > 0
        WaterRegion()
    elseif z >= -t
        LipidRegion()
    else
        CytosolRegion()
    end
end

"""
    V(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm, NMAX=1000)

Main interface that automatically determines regions and dispatches to appropriate method.
NMAX controls the convergence of infinite sums.
"""
function V(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm, NMAX=1000)
    charge = charge_region(h, t)
    eval = eval_region(z, t)
    V(z, charge, eval; ρ=ρ, t=t, h=h, NMAX=NMAX)
end
