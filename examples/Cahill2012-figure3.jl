# Cahill2012-figure3
# 
# Reproduces Figure 3 & 4 from:
# Cahill, K. (2012). Models of membrane electrostatics. Physical Review E, 85(5), 051921. 
# https://doi.org/10.1103/PhysRevE.85.051921

using MembraneElectrostatics
using Gnuplot
using ProgressMeter

#######
# Parameters
t = 5nm # membrane thickness - does this get passed through?

# "In the simulation, I let the potassium and chloride ions move according to a
# Metropolis algorithm within a box whose width and length were 50 nm and whose
# height was 10 nm. I took the potassium concentration to be 150 mM so as to
# allow for a 10 mM concentration of sodium ions. The box contained 2258 K+
# ions. The bottom of the box was covered by a uniform negative surface charge
# density whose total charge was − 143 |e| corresponding to 143 PSs at a molar
# density of  4%. I used 2115 Cl− ions to make the whole system neutral; these
# chloride ions played the role of the whole ensemble of anionic cell
# constituents."

# Create array of charges with XXX +1s and YYY -1s


# mini system, because of the O(N^2) cost of the explicit Coulomb pair sum #_#
charges=vcat(ones(90), -ones(85))
box=(10nm,10nm,10nm)

# midi system
# charges=vcat(ones(90*4), -ones(85*4))
# box=(20nm,20nm,10nm)

# Full system
# charges = vcat(ones(2258), -ones(2115))
# box = (50nm, 50nm, 10nm)

state = MembraneElectrostatics.MCState(charges, box)
# show(state)

# MembraneElectrostatics.calc_global_energy(state)
# show(state)

# The simulations consisted of eight separate runs in which 23_000 sweeps were
# allowed for thermalization. Four of the runs collected data for an additional
# 50_000 sweeps; the other four for an additional 9_000 sweeps.

SWEEPS=10_000

global ACCEPTED=0
@showprogress "MC sampling: " for i in 1:SWEEPS
    global ACCEPTED+=MembraneElectrostatics.mc_sweep!(state)
end
show(state)  # Show final state

println("Accepted $(ACCEPTED) Accepted ratio: $(ACCEPTED/(state.N*SWEEPS))")

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

# Figure 4

# Calculate average potential as function of z, for both ion types
# We'll sample z values from 0 to 10nm
z_values = range(0, 10e-9, length=100)

# Calculate potentials for K+ and Cl-
function test_charge_potential(state::MCState, z::Float64, test_charge::Float64)
    # Create temporary state with test charge at position (X,Y,z)
    temp_state = MCState(
        state.N + 1,
        hcat(state.positions, [state.L[1]/2, state.L[2]/2, z]),
        vcat(state.charges, test_charge),
        state.L,
        state.β,
        state.σ
    )
    
    # Calculate potential energy of test charge (last particle)
    V_no_self = calc_perion_energy(temp_state, temp_state.N, 
        CORRELATION=true, 
        ELECTROSTATIC=true, 
        SELF_INTERACTION=false) / test_charge

    V_with_self = calc_perion_energy(temp_state, temp_state.N,
        CORRELATION=true, 
        ELECTROSTATIC=true, 
        SELF_INTERACTION=true) / test_charge
    
    return V_no_self, V_with_self
end

# Calculate potentials at each z
V_K = zeros(length(z_values))
V_Cl = zeros(length(z_values))
V_no_self = zeros(length(z_values))

@showprogress "Calculating potentials: " for (i, z) in enumerate(z_values)
    V_no_self[i], V_K[i] = test_charge_potential(state, z, 1.0)
    _, V_Cl[i] = test_charge_potential(state, z, -1.0)
end

# Rescale to nm for plotting
z_nm = z_values ./ 1e-9

# Plot
@gp :- "set xlabel 'Distance from membrane (nm)'"
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set key right"

@gp z_nm V_K "with lines title 'K^+ potential' lw 2 lc rgb '#E41A1C'"
@gp :- z_nm V_Cl "with lines title 'Cl^- potential' dt 2 lw 2 lc rgb '#377EB8'"
@gp :- z_nm V_no_self "with lines title 'Without self-potential' dt 3 lw 2 lc rgb '#984EA3'"

Gnuplot.save("figure4.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure4.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")

# Hold my beer, I'm going to use printf

function debug_energy_components(S::MCState, i::Int)
    println("Ion $i analysis:")
    println("Position: ", S.positions[:,i])
    println("Charge: ", S.charges[i])
    
    E_corr = calc_perion_energy(S, i, CORRELATION=true, 
                               ELECTROSTATIC=false, SELF_INTERACTION=false)
    E_elec = calc_perion_energy(S, i, CORRELATION=false, 
                               ELECTROSTATIC=true, SELF_INTERACTION=false)
    E_self = calc_perion_energy(S, i, CORRELATION=false, 
                               ELECTROSTATIC=false, SELF_INTERACTION=true)
    
    println("Correlation energy: ", E_corr)
    println("Electrostatic energy: ", E_elec)
    println("Self-interaction: ", E_self)
end

for i in 1:state.N
    debug_energy_components(state, i)
end
