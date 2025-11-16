# Cahill2012-figure5
# 
# Reproduces Figure 5 from:
# Cahill, K. (2012). Models of membrane electrostatics. Physical Review E, 85(5), 051921. 
# https://doi.org/10.1103/PhysRevE.85.051921
#
# Figure 5 shows for a K+ ion:
# - V^w_w(z): self-interaction potential due to image charges
# - V_σ(z): potential due to surface charge field
# - Their sum: V^w_w(z) + V_σ(z)

using MembraneElectrostatics
using Gnuplot

#######
# Parameters
membrane = CAHILL_LIVER

# Surface charge density for 4% PS: 143 electrons over 50x50nm
σ_4pct = -143 * q / (50nm * 50nm)  # Negative charge

z_values = range(0, 10nm, length=1000)
z_nm = z_values ./ nm

# Calculate potentials
V_self = zeros(length(z_values))      # V^w_w(z) - self-interaction
V_sigma_4pct = zeros(length(z_values)) # V_σ(z) for 4% PS

for (i, z) in enumerate(z_values)
    # Self-interaction potential (image charges)
    V_self[i] = V(z, WaterRegion(), WaterRegion(), ρ=0.0, t=membrane.t, h=z, m=membrane, NMAX=20)
    
    # Surface charge potential: V_σ(z) = -σ|z|/(ε_w + ε_c)
    V_sigma_4pct[i] = -σ_4pct * abs(z) / (membrane.ϵ_w + membrane.ϵ_c)
end

# Since field is linear, 20% PS is just 5× the 4% value
V_sigma_20pct = 5 * V_sigma_4pct

# Sums
V_sum_4pct = V_self .+ V_sigma_4pct
V_sum_20pct = V_self .+ V_sigma_20pct

# Plot
@gp "set linetype 1 pi -1 pt 1 lc rgb '#E41A1C' dt solid # red"
@gp :- "set linetype 2 pi -1 pt 6 lc rgb '#377EB8' dt solid # blue"
@gp :- "set linetype 3 pi -1 pt 2 lc rgb '#4DAF4A' dt solid # green"
@gp :- "set linetype 4 pi -1 pt 4 lc rgb '#984EA3' dt solid # purple"
@gp :- "set linetype 5 pi -1 pt 1 lc rgb '#FF7F00' dt solid # orange"

@gp :- "set xlabel 'Distance from membrane z (nm)'"
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set key right"

# Plot curves for 4% PS
@gp :- z_nm V_sum_4pct "with lines title 'V^w_w(z) + V_σ(z) (4%)' lw 2 lc rgb '#E41A1C'"
@gp :- z_nm V_sigma_4pct "with lines title 'V_σ(z) (4%)' dt 2 lw 2 lc rgb '#377EB8'"

# Plot curves for 20% PS
@gp :- z_nm V_sum_20pct "with lines title 'V^w_w(z) + V_σ(z) (20%)' lw 2 lc rgb '#984EA3'"
@gp :- z_nm V_sigma_20pct "with lines title 'V_σ(z) (20%)' dt 2 lw 2 lc rgb '#FF7F00'"

# V_self is the same for both (doesn't depend on σ)
@gp :- z_nm V_self "with lines title 'V^w_w(z)' lw 2 lc rgb '#4DAF4A'"

@gp :- "set xrange [0:2]"
@gp :- "set yrange [0:0.07]"

Gnuplot.save("figure5.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure5.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")

println("Figure 5 saved to figure5.png and figure5.pdf")

