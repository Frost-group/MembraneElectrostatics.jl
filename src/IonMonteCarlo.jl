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
function MCState(charges::Vector, L; T::Float64=300.0)
    N = length(charges)
    # Initialize random positions within box
    positions = L .* rand(N, 3)' # scale all by scalar; or by individual (X,Y,Z)
    #  β = 1/(kB*T)
    β = 1.0 / (1.380649e-23 * T)  # kB in J/K

    # Calculate surface charge density
    #  From Cahill V, 143 electrons over 50x50nm surface = 4% PS
    σ = 143 * q / (50nm * 50nm)
    
    MCState(N, positions, charges, L, β, σ)
end

function calc_energy(S::MCState)
    E = 0.0
    
    # Screened Coulomb interaction between all pairs of ions
    for i in 1:S.N
        for j in (i+1):S.N # avoid double counting and self-interaction
            # Distance between i and j
            d=norm(S.positions[:,i] - S.positions[:,j])

            # TODO: Should include replicas in X&Y for PBCs!
            
            # Calculate potential energy between pair
            #  ASSUMES WATER BETWEEN ALL IONS
            #   What would Cahill do? (WWCD?)
            E += 1/ϵ_w * S.charges[i] * S.charges[j] * 1/d 
                # Uhm, are the units correct here? 
        end
    end

    # Electrostatic interaction between ions and membrane charge, see (27) in Cahill
    for i in 1:S.N
        E += S.positions[3,i] * S.charges[i] * S.σ / ϵ_w
    end

    # Now recurrance formulae; Eqn 9.
    for i in 1:S.N

    end

    return E
end

function mc_sweep!(S::MCState; δr=0.1nm)
    ACCEPTED = 0
    # One sweep = N attempted moves
    for i in 1:S.N # should we shuffle?

        # Store previous position and energy
        r_prev = S.positions[:,i]
        E_prev = calc_energy(S)
        
        # Trial move - randn displacement
        S.positions[:,i] += δr * randn(3)
        
        # Calculate new energy
        E_new = calc_energy(S)
        
        # Metropolis acceptance criterion
        ΔE = E_new - E_prev
        if ΔE < 0 || rand() < exp(-S.β * ΔE)
            ACCEPTED += 1 # Accept move
        else # Reject move & restore state
            S.positions[:,i] = r_prev
        end
        println("Energy: $E_new (ΔE = $ΔE)")
    end
end

