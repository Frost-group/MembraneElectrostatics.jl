module MembraneElectrostatics

using Gnuplot

# Physical constants
"Permittivity of free space, (C²N⁻¹m⁻²)."
const ε_0 = ϵ_0 = 8.85418682e-12

"Electron charge, (kgm²s⁻²)."
const eV = q = ElectronVolt = 1.602176634e-19

include("CahillRecursions.jl")
include("IonMonteCarlo.jl")

# Export constants
export ε_0, ϵ_0, eV, q, ElectronVolt, nm
# Export CahillRecursionsfunctions
export V_ww, V_wl, V_wc, V
# Export types
export WaterRegion, LipidRegion, CytosolRegion
export CahillMembrane, CAHILL_LIVER
# Export IonMonteCarlo functions
export MCState, calc_global_energy, calc_perion_energy


end 
