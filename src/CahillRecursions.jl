
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

# Helper functions
function summand9(n, z, ρ, t=5nm, h=1nm)
    p*p′^(n-1)/(√(ρ^2 + (z + 2n*t + h)^2))
end

function summand10(n, z, ρ, t=5nm, h=1nm)
    (p*p′)^n * (1/√(ρ^2 + (z - 2n*t - h)^2) - p′/√(ρ^2 + (z + 2*(n+1)*t + h)^2))
end

function summand11(n, z, ρ, t=5nm, h=1nm)
    (p*p′)^n / √(ρ^2 + (z - 2n*t - h)^2)
end

# Main potential functions
"""
    V_ww(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)

Calculate the potential in water according to Equation 9 from Cahill (2012).
"""
function V_ww(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)
    q/(4π*ϵ_w) * (
        1/√(ρ^2+(z-h)^2) 
        + p/√(ρ^2+(z+h)^2)
        - p′*(1-p^2) * sum(summand9(n,z,ρ,t,h) for n in 1:1000)
    )
end

"""
    V_wl(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)

Calculate the potential in the lipid bilayer according to Equation 10 from Cahill (2012).
"""
function V_wl(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)
    q/(4π*ϵ_wl) * sum(summand10(n,z,ρ,t,h) for n in 0:1000)
end

"""
    V_wc(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)

Calculate the potential in the cytosol according to Equation 11 from Cahill (2012).
"""
function V_wc(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)
    q*ϵ_l/(4π*ϵ_wl*ϵ_cl) * sum(summand11(n,z,ρ,t,h) for n in 0:1000)
end

"""
    V(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)

Calculate the potential at position z, automatically selecting the appropriate equation
based on the region (water, lipid, or cytosol).
"""
function V(z; ρ=sqrt(0.707nm^2+0.707nm^2), t=5nm, h=1nm)
    if z >= 0  # in water z>=0
        V_ww(z; ρ=ρ, t=t, h=h)
    elseif z > -t  # in lipid -t<z<0
        V_wl(z; ρ=ρ, t=t, h=h)
    else  # in cytosol, z<-t
        V_wc(z; ρ=ρ, t=t, h=h)
    end
end
