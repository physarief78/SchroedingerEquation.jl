using SchroedingerEquation
using LinearAlgebra
using CairoMakie

# Setup
basis = RealSpaceGrid1D(-50.0, 50.0, 1000)
# A barrier in the center (positive potential height 15)
barrier = SquareWellPotential(2.0, -15.0; x0=0.0) 
ham = build_hamiltonian(basis, barrier)

# Initial state: A moving Gaussian wavepacket starting on the left
x0 = -25.0
sigma = 2.0
p0 = 5.0 # E = p^2/2 = 12.5 < 15.0 (Barrier height) -> Tunneling regime
psi0_vals = [exp(-(x - x0)^2 / (2 * sigma^2)) * exp(im * p0 * x) for x in basis.x]
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)

# TDSE: Tunneling animation
t_total = 12.0
dt = 0.02
history = propagate_tdse(psi0, ham, t_total, dt; save_stride=5)

# Visualization
animate_evolution(history, barrier; filename="examples/results/03_tunneling.mp4", fps=30, title="Tunneling Effect")
println("Saved animation to examples/results/03_tunneling.mp4")
