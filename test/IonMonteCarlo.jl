

state = MembraneElectrostatics.MCState([1.0,1,1,-1,-1,-1], (10.0,10.0,10.0))

MembraneElectrostatics.calc_energy(state) == 0.0
