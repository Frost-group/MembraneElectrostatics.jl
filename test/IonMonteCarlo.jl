state = MembraneElectrostatics.MCState([1.0,1,1,-1,-1,-1], (10.0,10.0,10.0))

# Test that calc_perion_energy runs without error
energy = MembraneElectrostatics.calc_perion_energy(state, 1)
@test typeof(energy) == Float64

# Test that mc_sweep! runs without error
accepted = MembraneElectrostatics.mc_sweep!(state)
@test typeof(accepted) == Int64
@test accepted >= 0
