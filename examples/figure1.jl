using MembraneElectrostatics
using Gnuplot

# Parameters
t = 5nm
h = 1nm
Zs = collect(-(t+5nm):t/100:5nm)

# Calculate potentials
Vs = [V(z) for z in Zs]

# Plot
@gp "set xlabel 'Distance from membrane (nm)'"
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set arrow from 0,0 to 0,0.025 nohead lc rgb 'red'"
@gp :- "set arrow from -5E-9,0 to -5E-9,0.025 nohead lc rgb 'red'"
@gp :- Zs Vs "w l"


Gnuplot.save("figure1.png", term="png size 800,600 fontscale 0.8")
Gnuplot.save("figure1.pdf", term="pdf ")

