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
# Default liver membrane: ϵ_w=80*ϵ_0, ϵ_l=2*ϵ_0, ϵ_c=80*ϵ_0, t=5nm
membrane = CAHILL_LIVER

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
# ~3 mins (5nm cutoff) on my macpro with 23k sweeps
# 30 mins (Inf cutoff) ^_^
charges=vcat(+ones(90), -ones(85))
box=(10nm,10nm,10nm)

#charges=vcat(-ones(90), +ones(1))

# midi system ~ 20 mins (5 nm cutoff) on my macpro wiht 23k sweeps
#charges=vcat(ones(90*4), -ones(85*4))
#box=(20nm,20nm,10nm)

# Full system
#charges = vcat(ones(2258), -ones(2115))
#box = (50nm, 50nm, 10nm)

state = MembraneElectrostatics.MCState(charges, box)
# show(state)

# MembraneElectrostatics.calc_global_energy(state)
# show(state)

# "The simulations consisted of eight separate runs in which 23_000 sweeps were
# allowed for thermalization. Four of the runs collected data for an additional
# 50_000 sweeps; the other four for an additional 9_000 sweeps."

SWEEPS=23_000

global ACCEPTED=0
@showprogress "MC sampling: " for i in 1:SWEEPS
    global ACCEPTED+=MembraneElectrostatics.mc_sweep!(state, m=membrane)
end
#show(state)  # Show final state

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

# Plot histogram directly in Gnuplot
@gp :- "set style data histograms"
@gp :- "set style fill solid 0.5"
@gp :- "set boxwidth 0.1"
@gp :- "binwidth = 0.2"  # 10nm/50 bins
@gp :- "bin(x) = binwidth * floor(x/binwidth)"

# Separate z coordinates by charge
cation_z_nm = state.positions[3, state.charges .== 1.0] ./ 1e-9  # +1 charges
anion_z_nm = state.positions[3, state.charges .== -1.0] ./ 1e-9  # -1 charges

# Plot histograms with offset
@gp :- cation_z_nm "using (bin(\$1)):(1.0) smooth frequency with boxes title 'Cation density' lc rgb '#E41A1C' fs transparent solid 0.3"
@gp :- anion_z_nm "using (bin(\$1)+0.1):(1.0) smooth frequency with boxes title 'Anion density' lc rgb '#377EB8' fs transparent solid 0.3"

# Calculate means
cation_mean = length(cation_z_nm)/50
anion_mean = length(anion_z_nm)/50
# Add mean value lines
@gp :- "set arrow 1 from 0,$cation_mean to 10,$cation_mean nohead lc rgb '#E41A1C' lw 1"
@gp :- "set arrow 2 from 0,$anion_mean to 10,$anion_mean nohead lc rgb '#377EB8' lw 1"


cation_count = length(cation_z_nm)
anion_count = length(anion_z_nm)

# Add kernel density estimates
@gp :- "set samples 100"  # Increase smoothness of the kdensity
@gp :- cation_z_nm "using 1:(1.0/5) smooth kdensity bandwidth 0.2 with lines title 'Cation KDE' lw 2 lc rgb '#E41A1C'"
@gp :- anion_z_nm "using 1:(1.0/5) smooth kdensity bandwidth 0.2 with lines title 'Anion KDE' lw 2 lc rgb '#377EB8'"

Gnuplot.save("figure3.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure3.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")

# Figure 4

# Calculate average potential as function of z, for both ion types
# We'll sample z values from 0 to 10nm
z_values = range(0, 10e-9, length=100)

# Calculate potentials for K+ and Cl-
function test_charge_potential(state::MCState, z::Float64, test_charge::Float64; n=10, m=CAHILL_LIVER)
    x_points = range(0, state.L[1], length=n)
    y_points = range(0, state.L[2], length=n)
    
    V_full = 0.0        # correlation + self
    V_nocorr = 0.0      # no correlation + self
    V_noself = 0.0      # correlation + no self
    V_bare = 0.0        # no correlation + no self
    
    # Average over grid points
    @showprogress "Calculating potential at z=$(z/1e-9)nm: " for x in x_points, y in y_points
        temp_state = MCState(
            state.N + 1,
            hcat(state.positions, [x, y, z]),
            vcat(state.charges, test_charge),
            state.L,
            state.β,
            state.σ
        )
        
        # Calculate potential energy of test charge (last particle) in different ways
        V_full += calc_perion_energy(temp_state, temp_state.N, 
            CORRELATION=true, 
            ELECTROSTATIC=true, 
            SELF_INTERACTION=true, 
            m=m) / test_charge

        V_nocorr += calc_perion_energy(temp_state, temp_state.N,
            CORRELATION=false, 
            ELECTROSTATIC=true, 
            SELF_INTERACTION=true,
            m=m) / test_charge

        V_noself += calc_perion_energy(temp_state, temp_state.N,
            CORRELATION=true, 
            ELECTROSTATIC=true, 
            SELF_INTERACTION=false,
            m=m) / test_charge

        V_bare += calc_perion_energy(temp_state, temp_state.N,
            CORRELATION=false, 
            ELECTROSTATIC=true, 
            SELF_INTERACTION=false,
            m=m) / test_charge
    end
    
    # Return averaged potentials
    return V_full/(n*n), V_nocorr/(n*n), V_noself/(n*n), V_bare/(n*n)
end

# Calculate potentials at each z
V_K_full = zeros(length(z_values))
V_K_nocorr = zeros(length(z_values))
V_K_noself = zeros(length(z_values))
V_K_bare = zeros(length(z_values))

V_Cl_full = zeros(length(z_values))
V_Cl_nocorr = zeros(length(z_values))
V_Cl_noself = zeros(length(z_values))
V_Cl_bare = zeros(length(z_values))

@showprogress "Calculating potentials: " for (i, z) in enumerate(z_values)
    # K+ potentials
    V_K_full[i], V_K_nocorr[i], V_K_noself[i], V_K_bare[i] = test_charge_potential(state, z, 1.0, m=membrane)
    # Cl- potentials
    V_Cl_full[i], V_Cl_nocorr[i], V_Cl_noself[i], V_Cl_bare[i] = test_charge_potential(state, z, -1.0, m=membrane)
end

# Rescale to nm for plotting
z_nm = z_values ./ 1e-9

# Plot
@gp :- "set xlabel 'Distance from membrane (nm)'"
@gp :- "set ylabel 'Potential (V)'"
@gp :- "set key right"

# K+ potentials
@gp z_nm V_K_full "with lines title 'K^+ (full)' lw 2 lc rgb '#E41A1C'"
@gp :- z_nm V_K_nocorr "with lines title 'K^+ (no corr)' dt 2 lw 2 lc rgb '#E41A1C'"
@gp :- z_nm V_K_noself "with lines title 'K^+ (no self)' dt 3 lw 2 lc rgb '#E41A1C'"
@gp :- z_nm V_K_bare "with lines title 'K^+ (bare)' dt 4 lw 2 lc rgb '#E41A1C'"

# Cl- potentials
@gp :- z_nm V_Cl_full "with lines title 'Cl^- (full)' lw 2 lc rgb '#377EB8'"
@gp :- z_nm V_Cl_nocorr "with lines title 'Cl^- (no corr)' dt 2 lw 2 lc rgb '#377EB8'"
@gp :- z_nm V_Cl_noself "with lines title 'Cl^- (no self)' dt 3 lw 2 lc rgb '#377EB8'"
@gp :- z_nm V_Cl_bare "with lines title 'Cl^- (bare)' dt 4 lw 2 lc rgb '#377EB8'"

Gnuplot.save("figure4.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
Gnuplot.save("figure4.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")

# Hold my beer, I'm going to use printf

function debug_energy_components(S::MCState, i::Int; m=CAHILL_LIVER)
    println("Ion $i analysis:")
    println("Position: ", S.positions[:,i])
    println("Charge: ", S.charges[i])
    
    E_corr = calc_perion_energy(S, i, CORRELATION=true, 
                               ELECTROSTATIC=false, SELF_INTERACTION=false, m=m)
    E_elec = calc_perion_energy(S, i, CORRELATION=false, 
                               ELECTROSTATIC=true, SELF_INTERACTION=false, m=m)
    E_self = calc_perion_energy(S, i, CORRELATION=false, 
                               ELECTROSTATIC=false, SELF_INTERACTION=true, m=m)
    
    println("Correlation energy: ", E_corr)
    println("Electrostatic energy: ", E_elec)
    println("Self-interaction: ", E_self)
end

 for i in [1,2,state.N-1,state.N] # just print one +cation and one -anion
    debug_energy_components(state, i, m=membrane)
end
