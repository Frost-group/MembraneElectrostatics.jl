using LinearAlgebra

"""
Represents the state of a Monte Carlo simulation for an ionic system
"""
struct MCState
    N::Int64                    # Number of particles
    positions::Matrix{Float64}  # Nx3 matrix of particle positions
    charges::Vector{Float64}    # N-length vector of particle charges
    L::Tuple{Float64,Float64,Float64} # Lattice vectors / box size
    β::Float64                        # 1/(kB*T)
    σ::Float64                         # Surface charge density
end
"""
    MCState(charges::Vector{Float64}, L; T::Float64=310.15)

Initialize an (ion) Monte Carlo simulation state given charges and box size (scalar, or 3d tuple (X,Y,Z)).

Randomly distributes particles within the box.

Temperature T is in Kelvin (defaults to 310.15K = 37°C body temperature, matching Cahill's Fortran code).
"""
function MCState(charges::Vector, L; T::AbstractFloat=310.15)
    N = length(charges)
    # Initialize random positions within box
    positions = L .* rand(N, 3)' # scale random positions by scalar L; or by (X,Y,Z) tuple 
    #  β = 1/(kB*T)
    β = 1.0 / (8.6173303e-5 * T)  # kB in eV

    # Calculate surface charge density
    #  From Cahill V, 143 electrons over 50x50nm surface = 4% PS
    σ = 143 * q / (50nm * 50nm)
    # numerically about 1/100 ? Units C/m^2 ?
    
    MCState(N, positions, charges, L, β, σ)
end

# TODO: rewrite a function that just calculates 'atoms in molecules' sum for 'i' 
#  Nature of the loop suggests this should be fine. (Classical physics baby!)
function calc_global_energy(S::MCState)

#  WARNING: I think this is mostly WRONG!
    exception("calc_global_energy is not implemented")

end

#  Nature of the loop in the global energy function suggests this should be fine. 
#     (Classical physics baby!)
function calc_perion_energy(S::MCState, i::Int; cutoff=Inf, CORRELATION=true, ELECTROSTATIC=true, SELF_INTERACTION=true, m::CahillMembrane=CAHILL_LIVER)
    E = 0.0 # energy units of eV implicit everywhere
   
    qi=S.charges[i] # to multiply by all the calculation potentials V
    z=S.positions[3,i]

if CORRELATION 
    # Screened Coulomb interaction, mediated by the image charges induced in the membrane, between all pairs of ions, via Eqn 9
    # Hard-core potential parameters (Fortran kcl.f90 line 119: prevents K-Cl pairs from getting unphysically close)
    # KClsep = 0.236 nm minimum separation, hardcore ~ r^-12 repulsion
    KClsep = 0.236nm
    # Pre-factor for hard-core: in Fortran this is elec*KClsep^11/(12*4*pi*eps0*epsw)
    # We work in eV units, so: 
    hardCoreV_prefactor = (q * KClsep^11) / (12.0 * 4.0 * π * ε_0 * 80)    

    r_diff=Vector{Float64}(undef,3) # preallocate for loop
    for j in 1:S.N
        qj=S.charges[j]
        
        # Include central cell and nearest neighbors in x,y (9 boxes total: Fortran shift1, shift2 = -1,0,1)
        for dx in (-1,0,1), dy in (-1,0,1)
            # Skip self-interaction in central cell only (not in periodic images)
            if i==j && dx==0 && dy==0
                continue
            end
            
            # Offset positions by lattice vectors
            r_diff[1] = S.positions[1,i] - (S.positions[1,j] + dx*S.L[1])
            r_diff[2] = S.positions[2,i] - (S.positions[2,j] + dy*S.L[2])
            r_diff[3] = z - S.positions[3,j]

            ρ=sqrt(r_diff[1]^2 + r_diff[2]^2) # distance in plane
                            
            full_rsq = r_diff[1]^2 + r_diff[2]^2 + r_diff[3]^2
            if sqrt(full_rsq) > cutoff 
#               cutoff to reduce cost... but weird
#               artefacts due to the induced field? perhaps we should just restrict
#               to X,Y as before? 
                # tests suggest this DOOES introduce massive artefacts, needs to have smooth (erf?) cutoffs applied
                continue
            end
            
            # Safety check: skip unphysically close interactions (likely initialization artifact)
            # If ions are closer than KClsep/2, skip to prevent numerical blow-up in V
            if full_rsq < (KClsep / 2.0)^2
                continue
            end
            
            h=S.positions[3,j] # vertical position of jth ion
            
            #   What would Cahill do? (WWCD?)
            #   "I went toward you, endlessly toward the light"
 
            # _Potential_ from image charges generated in membrane by the charge at r_diff
            V_elec = V(z,WaterRegion(),WaterRegion(), ρ=ρ, t=m.t, h=h, m=m, NMAX=10)
            # This is definitely in Volts, as all the figures reproduce. 
            
            # Add hard-core repulsion for opposite-charge pairs (K-Cl)
            # Fortran kcl.f90 lines 194, 261: EKCl = EKCl + elec*(Volt + hardCoreV/rsq**6)
            if qi * qj < 0  # Opposite charges
                # Safety check: prevent numerical blow-up if ions are unphysically close
                # Cap rsq at KClsep^2 to ensure hard-core is large but finite
                # This prevents blow-up while still giving strong repulsion for close ions
                rsq_min = KClsep^2
                rsq = max(full_rsq, rsq_min)
                # Hard-core contribution: V_hc = prefactor / r^12 (in Volts)
                V_hardcore = hardCoreV_prefactor / (rsq^6)
                E += qi * qj * (V_elec + V_hardcore) # q is (-1,+1) , so unitwise still eV
            else
                E += qi * qj * V_elec
            end
        end
    end
end

# Electrostatic interaction between ions and membrane charge, see (27) in Cahill
# constant E-field -> potential V = z * Ef = z * σ / (ϵ_w+ϵ_c)
# units of eV presumed ? Is that the correct dielectric constant?
if ELECTROSTATIC
    E += qi * (z * S.σ / (m.ϵ_w+m.ϵ_c) ) 
end

# recurrance formulae for ion self-interaction with slab dielectrics (membrane); Eqn 9.
# now detects ρ≈0 and runs the dedicated (no infinite 1/r term) code
if SELF_INTERACTION
    E+= qi * qi * V(z,WaterRegion(),WaterRegion(), ρ=0.0, t=m.t, h=z, m=m, NMAX=20)
end

    return E
end

function mc_sweep!(S::MCState; δr=0.5nm, GLOBAL_ENERGY=false, m::CahillMembrane=CAHILL_LIVER)
    ACCEPTED = 0 # One sweep = N attempted moves
    height = S.L[3]

    for i in 1:S.N # should we shuffle?
        r = @view S.positions[:, i]
        
        # Store previous position and energy
        # tuple is stack-allocated, no heap allocation
        r_prev = (r[1], r[2], r[3])
        E_prev = GLOBAL_ENERGY ? calc_global_energy(S) : calc_perion_energy(S, i, m=m)
        
        # Trial move - randn displacement
        r[1] += δr * randn()
        r[2] += δr * randn()
        r[3] += δr * randn()
        
        # Now follow the Fortran code to use reflection in Z; 
        # stop ions from penetrating membrane / leaving box
        # Fortran code (kcl.f90 lines 237-242):
        #       if z > height: z = 2*height - z 
        #       if z < 0: z = -z
        r[3] < 0 && (r[3] = -r[3])
        r[3] > height && (r[3] = 2height - r[3])

        # PBCs in X and Y; project back into box
        r[1] = mod(r[1], S.L[1])
        r[2] = mod(r[2], S.L[2])
        
        # Calculate new energy
        E_new = GLOBAL_ENERGY ? calc_global_energy(S) : calc_perion_energy(S, i, m=m)
        
        # Metropolis acceptance criterion
        ΔE = E_new - E_prev
        if ΔE < 0 || rand() < exp(-S.β * ΔE)
            ACCEPTED += 1 # Accept move
        else # Reject move & restore state
            r[1], r[2], r[3] = r_prev
        end
#        println("Energy: $E_new (ΔE = $ΔE)")
    end

    return ACCEPTED
end

