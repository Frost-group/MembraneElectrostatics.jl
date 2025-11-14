"""
Cahill2012-Figure78.jl
Transfer Matrix Examples: Phospholipid Bilayer with Head Groups

Reproduces Cahill's 5-region model (3 internal layers, phosopholipds with crazy
high permittivity), in comparison to the 'Liver' membrane model.  
"""

using MembraneElectrostatics
using Gnuplot
using Printf
using Statistics

# Load transfer matrix module
include("../src/TransferMatrix.jl")

# Physical constants
const nm = 1e-9
const ϵ_0 = 8.85418682e-12
const q = 1.602176634e-19

# Calculate 4% PS surface charge density
# Cahill: 143 electrons over 50×50 nm surface at 4% PS
σ_4percent = -143 * q / (50nm * 50nm)  # C/m², negative for PS



println("=" ^ 70)
println("Phospholipid Bilayer with Head Groups: Transfer Matrix Method")
println("=" ^ 70)

# System 1: Bare lipid slab (1 layer, 3 regions total)
# Cahill's 'Liver' membrane model: single slab
println("\nSystem 1: Bare lipid bilayer (5 nm)")
println("  Water (ϵ=80ϵ₀) | Lipid (ϵ=2ϵ₀, 5nm) | Cytosol (ϵ=80ϵ₀)")

ϵ_w = 80 * ϵ_0
ϵ_lipid = 2 * ϵ_0
ϵ_c = 80 * ϵ_0

lipid_only = [TransferMatrix.Layer(ϵ_lipid, 5 * nm)]
system_1layer = TransferMatrix.LayeredSystem(ϵ_w, ϵ_c, lipid_only)

# System 2: Lipid bilayer with head groups (3 layers, 5 regions total)
# Cahill Section "Three Dielectric Layers"
# Based on Stern & Feller (2003), Nymeyer & Zhou (2008), Baker (2011)
#  but actually Cahill 'simplifies' it to what he fancies. I suppose now the code is automated we could go back and do the full thing? 
println("\nSystem 2: Phospholipid bilayer with head groups (5 nm total), Cahill bettween Eqn 55 and 56")
println("  Water (ϵ=80ϵ₀) | Head₁ (ϵ=195ϵ₀, 1nm) | Lipid (ϵ=2ϵ₀, 3nm) | Head₂ (ϵ=195ϵ₀, 1nm) | Cytosol (ϵ=80ϵ₀)")

ϵ_head = 195 * ϵ_0  # High permittivity phosphate head groups
layers_3 = [
    TransferMatrix.Layer(ϵ_head, 1 * nm),   # Outer head group
    TransferMatrix.Layer(ϵ_lipid, 3 * nm),  # Lipid core
    TransferMatrix.Layer(ϵ_head, 1 * nm)    # Inner head group
]

system_3layer = TransferMatrix.LayeredSystem(ϵ_w, ϵ_c, layers_3)

#######
# Fig1b: Potential vs. height, via transmfer matrix and normal recursions 
#######

zs = range(0.0, 5.0 * nm, length=50) # length = number of points in curve 

V_transfer_matrix = [TransferMatrix.V_transfer_matrix(1nm, z, 1nm, system_1layer, q=q) for z in zs]
V_recursions = [V(z, ρ=1nm, h=1nm, m=CAHILL_LIVER, NMAX=1000) for z in zs]

@gp :- "set xlabel 'Z - height above membrane (nm)'"
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set key top right"
@gp :- "set yrange [0:0.03]" # a la Cahill 

@gp :- zs./nm V_transfer_matrix "w l lw 2 lc rgb '#984EA3' title 'Transfer Matrix'"
@gp :- zs./nm V_recursions "w l lw 3 lc rgb '#377EB8' title 'Recursions'"

Gnuplot.save("figure1b.png", 
             term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure1b.pdf",
             term="pdfcairo size 4in,3in enhanced font 'Helvetica,10'")
println("  Saved: figure1b.png")
println("  Saved: figure1b.pdf")

#######
# Fig7: Potential vs. height, with a charge at h=1nm, ρ=1nm for both systems
#######

# charge at h=1nm, ρ=1nm, evaluate from z=0 to z=5nm
h = 1.0 * nm  # Charge height above membrane (in water)
ρ = 1.0 * nm  # Radial offset

zs = range(0.0, 5.0 * nm, length=50) # length = number of points in curve 

println("\nComputing potentials...")
println("  Charge: q = |e| at (ρ,h) = ($(ρ/nm), $(h/nm))")
println("  Evaluating potential from z=0 to z=$(zs[end]/nm)nm")

# Compute potentials (neutral and charged membrane)
V_1layer = Float64[]
V_3layer = Float64[]
V_1layer_PS = Float64[]
V_3layer_PS = Float64[]

print("  Progress: ")
for (i, z) in enumerate(zs)
    # 1-layer system (bare lipid)
    V1 = TransferMatrix.V_transfer_matrix(ρ, z, h, system_1layer, q=q)
    push!(V_1layer, V1)
    
    # 3-layer system (with head groups)
    V3 = TransferMatrix.V_transfer_matrix(ρ, z, h, system_3layer, q=q)
    push!(V_3layer, V3)
    
    # Add 4% PS surface charge contribution
    # Cahill Eq. 27: V_PS(z) = -z · σ / (ϵ_w + ϵ_c)
    V_PS = -z * σ_4percent / (ϵ_w + ϵ_c)
    push!(V_1layer_PS, V1 + V_PS)
    push!(V_3layer_PS, V3 + V_PS)
    
    if i % 10 == 0 # this was a bit slow, but faster now with the cutoffs and early returns 
        print("$(i)/$(length(zs))... ")
    end
end
println("Done!")

# Convert to plotting units
z_nm = zs ./ nm
V_1layer_V = V_1layer  # Already in Volts
V_3layer_V = V_3layer

#######
# Plot: Cahill Figure 7
#######

println("\nGenerating plots...")

# Setup Gnuplot color scheme (Cahill's colors)
@gp "set linetype 1 pi -1 pt 1 lc rgb '#E41A1C' dt solid # red"
@gp :- "set linetype 2 pi -1 pt 6 lc rgb '#377EB8' dt solid # blue"
@gp :- "set linetype 3 pi -1 pt 2 lc rgb '#4DAF4A' dt solid # green"
@gp :- "set linetype 4 pi -1 pt 4 lc rgb '#984EA3' dt solid # magenta"
@gp :- "set linetype 5 pi -1 pt 1 lc rgb '#FF7F00' dt solid # orange"

# Main plot
@gp :- "set xlabel 'Z - height above membrane (nm)'"
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set key top right"

# Mark membrane boundaries
@gp :- "set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'black' lw 1 dt 2"
@gp :- "set label 'Membrane' at 0.2, graph 0.95 left"

# Plot neutral membrane curves
@gp :- z_nm V_1layer_V "w l lw 2 lc rgb '#984EA3' dt 2 title 'Bare lipid (neutral)'"
@gp :- z_nm V_3layer_V "w l lw 3 lc rgb '#377EB8' dt 2 title 'With head groups (neutral)'"

# Plot charged membrane curves (4% PS)
@gp :- z_nm V_1layer_PS "w l lw 2 lc rgb '#E41A1C' title 'Bare lipid (4% PS)'"
@gp :- z_nm V_3layer_PS "w l lw 3 lc rgb '#FF7F00' title 'With head groups (4% PS)'"

# Add reference: direct Coulomb for comparison
r_direct = sqrt.(ρ^2 .+ (zs .- h).^2)
V_coulomb = q ./ (4π * ϵ_w .* r_direct)
@gp :- z_nm V_coulomb "w l lw 1 lc rgb '#4DAF4A' dt 3 title 'Pure Coulomb, ϵ_w, (no membrane effects)'"

@gp :- "set yrange [0:0.03]" # a la Cahill 

# Save plots
Gnuplot.save("figure7.png", 
             term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure7.pdf",
             term="pdfcairo size 4in,3in enhanced font 'Helvetica,10'")

println("  Saved: figure7.png")
println("  Saved: figure7.pdf")

#######
# Ion-interaction energy with membranes, with and without 4% PS surface charge
#######

println("\n" * "=" ^ 70)
println("Self-Interaction Energy: Ion approaching membrane")
println("=" ^ 70)


# need to avoid infinity at z=0
zs = range(0.1nm, 5.0nm, length=100) # length = number of points in curve 

# Arrays for neutral membrane
V_self_1layer = Float64[]
V_self_3layer = Float64[]

# Arrays for charged membrane (4% PS)
V_self_1layer_PS = Float64[]
V_self_3layer_PS = Float64[]



println("\nComputing self-interaction potentials...")
println("  Surface charge density (4% PS): σ = $(σ_4percent) C/m²")

for z in zs
    # Image charge contributions (neutral membrane)
    V1 = TransferMatrix.V_transfer_matrix_self_interaction(z, system_1layer, q=q)
    V3 = TransferMatrix.V_transfer_matrix_self_interaction(z, system_3layer, q=q)
    
    # Neutral membrane
    push!(V_self_1layer, V1)
    push!(V_self_3layer, V3)
    
    # Charged membrane: add uniform field from PS surface charge
    # Cahill Eq. 27: V_PS(z) = -z · σ / (ϵ_w + ϵ_c)
    V_PS = -z * σ_4percent / (ϵ_w + ϵ_c)
    
    push!(V_self_1layer_PS, V1 + V_PS)
    push!(V_self_3layer_PS, V3 + V_PS)
end

z_self_nm = zs ./ nm

# Plot self-interaction potentials
@gp "set xlabel 'Z - height above membrane (nm)'"
@gp :- "set ylabel 'Self-interaction potential (V)'"
@gp :- "set title 'Fig8: Cation Potential Barrier'"
@gp :- "set key top right"

@gp :- "set yrange [-0.03:0.03]" # following Cahill, he let's the 'no PSs' curve go off the top

# Neutral membrane curves
@gp :- z_self_nm V_self_1layer "w l lw 2 lc rgb '#984EA3' dt 2 title 'Bare lipid (neutral)'"
@gp :- z_self_nm V_self_3layer "w l lw 3 lc rgb '#377EB8' dt 2 title 'With head groups (neutral)'"

# Charged membrane curves (4% PS)
@gp :- z_self_nm V_self_1layer_PS "w l lw 2 lc rgb '#E41A1C' title 'Bare lipid (4% PS)'"
@gp :- z_self_nm V_self_3layer_PS "w l lw 3 lc rgb '#FF7F00' title 'With head groups (4% PS)'"

Gnuplot.save("figure8-cation.png",
             term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure8-cation.pdf",
             term="pdfcairo size 4in,3in enhanced font 'Helvetica,10'")

println("\n  Saved: figure8-cation.png")
println("  Saved: figure8-cation.pdf")

