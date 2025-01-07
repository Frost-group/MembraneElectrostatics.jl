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
    MCState(charges::Vector{Float64}, L; T::Float64=300.0)

Initialize an (ion) Monte Carlo simulation state given charges and box size (scalar, or 3d tuple (X,Y,Z)).

Randomly distributes particles within the box.

Temperature T is in Kelvin (defaults to 300K).
"""
function MCState(charges::Vector, L; T::AbstractFloat=300.0)
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
function calc_perion_energy(S::MCState, i::Int; cutoff=2nm, CORRELATION=true, ELECTROSTATIC=true, SELF_INTERACTION=true)
    E = 0.0 # energy units of eV implicit everywhere
   
    qi=S.charges[i] # to multiply by all the calculation potentials V
    z=S.positions[3,i]

if CORRELATION 
    # Screened Coulomb interaction, mediated by the image charges induced in the membrane, between all pairs of ions, via Eqn 9
    r_diff=Vector{Float64}(undef,3) # preallocate for loop
    for j in 1:S.N
        if i==j # avoid self-interaction
            continue
        end 

        qj=S.charges[j]

        # Include central cell and nearest neighbors in x,y; naive approach therefore gives a 9x slowdown
        for dx in (-1,0,1), dy in (-1,0,1)
            # Offset positions by lattice vectors
            r_diff[1] = S.positions[1,i] - (S.positions[1,j] + dx*S.L[1])
            r_diff[2] = S.positions[2,i] - (S.positions[2,j] + dy*S.L[2])

            ρ=sqrt(r_diff[1]^2 + r_diff[2]^2) # distance in plane
            #  implement a cutoff...
            if ρ > cutoff
                continue
            end

            h=S.positions[3,j] # vertical position of jth ion
            #   What would Cahill do? (WWCD?)
            #   "I went toward you, endlessly toward the light"
            # Potential from image charges generated in membrane by the charge at r_diff
            E+=qi * qj * V(z,WaterRegion(),WaterRegion(), ρ=ρ, t=5nm, h=h, NMAX=100) 
                # Uhm, are the units correct here? 
        end
    end
end

# Electrostatic interaction between ions and membrane charge, see (27) in Cahill
# constant E-field -> potential V = z * Ef = z * σ / ϵ_w
# units of eV presumed ? Is that the correct dielectric constant?
if ELECTROSTATIC
    E += qi * (z * S.σ / ϵ_w ) 
end

# recurrance formulae for ion self-interaction with slab dielectrics (membrane); Eqn 9.
# now detects ρ≈0 and runs the dedicated (no infinite 1/r term) code
if SELF_INTERACTION
    E+=qi * qi * V(z,WaterRegion(),WaterRegion(), ρ=0.0, t=5nm, h=z, NMAX=100)
end

    return E
end

function mc_sweep!(S::MCState; δr=0.5nm, GLOBAL_ENERGY=false)
    ACCEPTED = 0 # One sweep = N attempted moves

    for i in 1:S.N # should we shuffle?

        # Store previous position and energy
        r_prev = S.positions[:,i]
        E_prev = GLOBAL_ENERGY ? calc_global_energy(S) : calc_perion_energy(S,i)
        
        # Trial move - randn displacement
        S.positions[:,i] += δr * randn(3)
        
        # Hard wall in Z; stop ions from penetrating membrane / leaving 10 nm box (magic number, FIXME)
        if S.positions[3,i] < 0.0 # stop ions from penetrating membrane
            S.positions[3,i] = 0.0
        end
        if S.positions[3,i] > 10.0nm
            S.positions[3,i] = 10.0nm # stop them escaping
        end

        # PBCs in X and Y; project back into box
        S.positions[1,i] = mod(S.positions[1,i], S.L[1])
        S.positions[2,i] = mod(S.positions[2,i], S.L[2])
        
        # Calculate new energy
        E_new = GLOBAL_ENERGY ? calc_global_energy(S) : calc_perion_energy(S,i)
        
        # Metropolis acceptance criterion
        ΔE = E_new - E_prev
        if ΔE < 0 || rand() < exp(-S.β * ΔE)
            ACCEPTED += 1 # Accept move
        else # Reject move & restore state
            S.positions[:,i] = r_prev
        end
#        println("Energy: $E_new (ΔE = $ΔE)")
    end

    return ACCEPTED
end

