# Fig1Potentials_Bacteria.jl
#
# Cahill, K. (2012). Models of membrane electrostatics. Physical Review E, 85(5), 051921.
# https://doi.org/10.1103/PhysRevE.85.051921

using MembraneElectrostatics
using Gnuplot

# Parameters
membranes = [CahillMembrane(t=5nm), CahillMembrane(t=10nm)]

h = 1nm

maxthickness = maximum([m.t for m in membranes])
Zs = collect(-(maxthickness+5nm):maxthickness/100:5nm)

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
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set arrow from 0,0 to 0,0.025 nohead lc rgb 'red'"
for membrane in membranes
    @gp :- "set arrow from -$(membrane.t),0 to -$(membrane.t),0.025 nohead lc rgb 'red'"
end
@gp :- "set xrange [-$(maxthickness+5nm):$(5nm)]"

# Calculate potentials
for membrane in membranes
    for h in [1nm]
        VCs = [V(z, h=h, m=membrane, NMAX=10) for z in Zs]
        @gp :- Zs VCs "w l title 'membrane thickness=$(membrane.t/nm) nm,  h=$(h/nm) nm'"
    end
end


Gnuplot.save("Fig1Potentials_Bacteria.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("Fig1Potentials_Bacteria.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")

