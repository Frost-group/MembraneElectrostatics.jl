# Cahill2012-figure3
#
# Cahill, K. (2012). Models of membrane electrostatics. Physical Review E, 85(5), 051921. 
# https://doi.org/10.1103/PhysRevE.85.051921

using MembraneElectrostatics
using Gnuplot
#using ProgressMeter # Pkg is a bit borked on my laptop, so I can't install this

# "In the simulation, I let the potassium and chloride ions move according to a
# Metropolis algorithm within a box whose width and length were 50 nm and whose
# height was 10 nm. I took the potassium concentration to be 150 mM so as to
# allow for a 10 mM concentration of sodium ions. The box contained 2258 K+
# ions. The bottom of the box was covered by a uniform negative surface charge
# density whose total charge was − 143 |e| corresponding to 143 PSs at a molar
# density of  4%. I used 2115 Cl− ions to make the whole system neutral; these
# chloride ions played the role of the whole ensemble of anionic cell
# constituents."

# Create array of charges with 2258 +1s and 2115 -1s
# charges = vcat(ones(2258), -ones(2115))
# Box 50nm x 50nm x 10nm
# box = (50nm, 50nm, 10nm)

#######
# Parameters
t = 5nm # membrane thickness - does this get passed through?
# mini system, because of the O(N^2) cost of the explicit Coulomb pair sum #_#
charges=vcat(ones(90), -ones(85))
box=(10nm,10nm,10nm)

state = MembraneElectrostatics.MCState(charges, box)
show(state)

MembraneElectrostatics.calc_global_energy(state)
show(state)

#@showprogress "MC sampling: " 
for i in 1:10_000
    print(".")
    MembraneElectrostatics.mc_sweep!(state)
end
show(state)  # Show final state

# Check convergence ¯\_(ツ)_/¯


# rescued from https://github.com/jarvist/gnuplot-snippets/blob/master/gnuplot-render.gpt 
@gp "set linetype 1 pi -1 pt 1 lc rgb '#E41A1C' dt solid # red"
@gp :- "set linetype 2 pi -1 pt 6 lc rgb '#377EB8' dt solid # blue"
@gp :- "set linetype 3 pi -1 pt 2 lc rgb '#4DAF4A' dt solid # green"
@gp :- "set linetype 4 pi -1 pt 4 lc rgb '#984EA3' dt solid # purple"
@gp :- "set linetype 5 pi -1 pt 1 lc rgb '#FF7F00' dt solid # orange"
@gp :- "set linetype 6 pi -1 pt 6 lc rgb '#FFFF33' dt solid # yellow"
@gp :- "set linetype 7 pi -1 pt 2 lc rgb '#A65628' dt solid # brown"
@gp :- "set linetype 8 pi -1 pt 4 lc rgb '#F781BF' dt solid # pink"


# Plot
@gp :- "set xlabel 'Distance from membrane (nm)'"
@gp :- "set ylabel 'Density'"
@gp :- "set arrow from 0,0 to 0,0.025 nohead lc rgb 'red'"
@gp :- "set arrow from -5E-9,0 to -5E-9,0.025 nohead lc rgb 'red'"

# Plot histogram directly in Gnuplot
@gp :- "set style data histograms"
@gp :- "set style fill solid 0.5"
@gp :- "set boxwidth 0.2"
@gp :- "binwidth = 0.2"  # 10nm/50 bins
@gp :- "bin(x) = binwidth * floor(x/binwidth)"

# Separate z coordinates by charge
cation_z_nm = state.positions[3, state.charges .== 1.0] ./ 1e-9  # +1 charges
anion_z_nm = state.positions[3, state.charges .== -1.0] ./ 1e-9  # -1 charges

# Plot both populations separately
@gp :- cation_z_nm "using (bin(\$1)):(1.0) smooth frequency with boxes title 'Cation density' lc rgb '#E41A1C'"
@gp :- anion_z_nm "using (bin(\$1)):(1.0) smooth frequency with boxes title 'Anion density' lc rgb '#377EB8'"

Gnuplot.save("figure3.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure3.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")

