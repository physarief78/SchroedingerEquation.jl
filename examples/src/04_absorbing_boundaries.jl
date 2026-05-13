using SchroedingerEquation
using LinearAlgebra
using CairoMakie

# Setup
basis = RealSpaceGrid1D(-20.0, 20.0, 1000) # Higher resolution for smoother CAP
# Use a slightly larger and stronger CAP for "0% reflections"
cap = AbsorbingBoundary(8.0, 2.0)
ham = build_hamiltonian(basis, CustomPotential(x -> 0.0), cap)

# Wavepacket
x0 = -5.0         # Start closer to center to give it space
p0 = 3.0          # Momentum
direction = -1.0  # Moving to the left towards the CAP

psi0_vals = exp.(-(basis.x .- x0).^2 ./ 2.0) .* exp.(im * p0 * direction .* basis.x)
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)

# Track norm
obs = (norm = wf -> sum(abs2, wf.psi) * wf.basis.dx,)
history, obs_data = propagate_tdse(psi0, ham, 20.0, 0.05; observables=obs)

# Visualization
fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], title="Survival Probability (CAP Absorption)", xlabel="Time Step", ylabel="Norm")
lines!(ax, obs_data.norm, color=:red, linewidth=2.5)
save("examples/results/04_absorbing_boundaries.png", fig)

# Animate without passing the boundary object as a potential
animate_evolution(history; filename="examples/results/04_absorbing.mp4", fps=30, title="Absorption Dynamics")
println("Saved results to examples/results/04_absorbing_boundaries.png and .mp4")
