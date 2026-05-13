using SchroedingerEquation
using LinearAlgebra
using GLMakie

# Setup
basis = RealSpaceGrid1D(-20.0, 20.0, 1000)
# Tunneling barrier
potential = SquareWellPotential(2.0, -10.0; x0=0.0) 
ham = build_hamiltonian(basis, potential)

# Wavepacket
x0 = -10.0
sigma = 2.0
p0 = 3.0
psi0_vals = [exp(-(x - x0)^2 / (2 * sigma^2)) * exp(im * p0 * x) for x in basis.x]
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)

println("Starting live interactive simulation...")
live_evolution(psi0, ham, 30.0, 0.05; stride=1)
println("Live simulation closed.")
 