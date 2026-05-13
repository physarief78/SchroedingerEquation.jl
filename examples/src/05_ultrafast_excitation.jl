using SchroedingerEquation
using LinearAlgebra
using SparseArrays
using CairoMakie

# Setup
basis = RealSpaceGrid1D(-30.0, 30.0, 1000)
static_pot = HarmonicPotential(1.0)
boundary = AbsorbingBoundary(5.0)
# Laser pulse
pulse = CustomTimeDependentPotential((x, t) -> -x * 0.1 * exp(-((t - 5.0)/1.0)^2) * cos(1.0 * t))
# Ground state
ham_static = build_hamiltonian(basis, static_pot, boundary)
energies, wavefunctions = solve_tise(ham_static, 1)
psi0 = wavefunctions[1]

# Dynamic propagation
# Note: We pass the full static Hamiltonian matrix here. 
# The solver will add the time-dependent 'pulse' to the diagonal of this matrix.
t_total = 12.0
dt = 0.05
# Track normalized energy expectation value
obs = (E = wf -> real(expectation_value(wf, ham_static.H_matrix)) / (sum(abs2, wf.psi) * wf.basis.dx),)
history, obs_data = propagate_tdse(psi0, basis, ham_static.H_matrix, pulse, t_total, dt; observables=obs)

# Visualization
fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], title="Energy Expectation Value <E(t)>", xlabel="Time Step", ylabel="Energy")
lines!(ax, real.(obs_data.E), color=:purple, linewidth=2.5)
save("examples/results/05_excitation.png", fig)
println("Saved results to examples/results/05_excitation.png")
