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
    E = 0.0 # eV implicit everywhere
    
    # Screened Coulomb interaction between all pairs of ions
    E_C=0.0
    r_diff=Vector{Float64}(undef,3) # preallocate for loop
    for i in 1:S.N
        for j in (i+1):S.N # avoid double counting and self-interaction
            # Distance between i and j
            r_diff .= view(S.positions,:,i) - view(S.positions,:,j) # views to avoid slices
            d=norm(r_diff)

            # TODO: Should include replicas in X&Y for PBCs!
            
            # Calculate potential energy between pair
            #  ASSUMES WATER BETWEEN ALL IONS
            #   What would Cahill do? (WWCD?)
            E_C += S.charges[i] * S.charges[j] / d 
                # Uhm, are the units correct here? 
        end
    end
    E *= E_C* q/(4π*ϵ_w) 

    # Electrostatic interaction between ions and membrane charge, see (27) in Cahill
    for i in 1:S.N
        E += S.positions[3,i] * S.charges[i] * S.σ / ϵ_w
        # Constant E-field -> potential proport to z, units of eV
    end

    # Now recurrance formulae; Eqn 9.
    for i in 1:S.N
        # which one to use?
        #  what about other charges interacting with the images charges? 
        #  Or would that be double counting?
        # function V(z, charge::WaterRegion, eval::WaterRegion; ρ, t, h, NMAX=1000)
        z=S.positions[3,i]
        E+=V(z,WaterRegion(),WaterRegion(), ρ=0.0, t=5nm, h=z, NMAX=100)
    end

    return E
end

#  Nature of the loop in the global energy function suggests this should be fine. 
#     (Classical physics baby!)
function calc_perion_energy(S::MCState, i::Int)
    E = 0.0 # eV implicit everywhere
   
    qi=S.charges[i] # to multiply by all the calculation potentials V
    
    # Screened Coulomb interaction between all pairs of ions, via Eqn 9
    E_C=0.0
    r_diff=Vector{Float64}(undef,3) # preallocate for loop
    for j in 1:S.N
        if i==j # avoid self-interaction
            continue
        end 

        # Distance between i and j
        r_diff .= view(S.positions,:,i) - view(S.positions,:,j) # views to avoid slices
        z=S.positions[3,i]
        ρ=sqrt(r_diff[1]^2 + r_diff[2]^2) 

        # TODO: Should include replicas in X&Y for PBCs!
        
        # Simple electrostatics
        #E_C += S.charges[i] * S.charges[j] / d

        #   What would Cahill do? (WWCD?)
        #   "I went toward you, endlessly toward the light"
        E_C+=qi * V(z,WaterRegion(),WaterRegion(), ρ=ρ, t=5nm, h=z, NMAX=100) 
            # Uhm, are the units correct here? 
    end
#    E *= E_C* q/(4π*ϵ_w)         #  ASSUMES WATER BETWEEN ALL IONS

    # Electrostatic interaction between ions and membrane charge, see (27) in Cahill
    E += S.positions[3,i] * qi * S.σ / ϵ_w 

# recurrance formulae for self-interaction with slab dielectrics (membrane); Eqn 9.
    z=S.positions[3,i]
    # currently eval by taking ρ to very large... I think this is the same as Eqn. 35
    #   FIXME: Actually implement the more simple Eqn. 35 (And maybe check understanding at same time) 
    E+=qi * V(z,WaterRegion(),WaterRegion(), ρ=100.0, t=5nm, h=z, NMAX=100)

end

function mc_sweep!(S::MCState; δr=0.5nm, GLOBAL_ENERGY=false)
    ACCEPTED = 0
    # One sweep = N attempted moves
    for i in 1:S.N # should we shuffle?

        # Store previous position and energy
        r_prev = S.positions[:,i]
        E_prev = GLOBAL_ENERGY ? calc_global_energy(S) : calc_perion_energy(S,i)
        
        # Trial move - randn displacement
        S.positions[:,i] += δr * randn(3)
        
        if S.positions[3,i] < 0.0 # stop ions from penetrating membrane
            S.positions[3,i] = 0.0
        end
        if S.positions[3,i] > 10.0nm
            S.positions[3,i] = 10.0nm # stop them escaping
        end
        
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

