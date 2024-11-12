module MembraneElectrostatics

using Gnuplot

# Physical constants
"Permittivity of free space, (C²N⁻¹m⁻²)."
const ε_0 = ϵ_0 = 8.85418682e-12

"Electron charge, (kgm²s⁻²)."
const eV = q = ElectronVolt = 1.602176634e-19


# Export constants
export ε_0, ϵ_0, eV, q, ElectronVolt

# Material constants
export ϵ_w, ϵ_l, ϵ_c, nm

# Export functions
export V_ww, V_wl, V_wc, V

include("CahillRecursions.jl")

end 