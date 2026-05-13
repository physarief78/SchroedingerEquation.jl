using SchroedingerEquation
using LinearAlgebra
using CairoMakie

# Setup
basis = RealSpaceGrid1D(-20.0, 20.0, 500)
potential = HarmonicPotential(1.0)
ham = build_hamiltonian(basis, potential)

# Initial state (coherent state/displaced gaussian)
psi0_vals = exp.(-(basis.x .+ 5.0).^2 ./ 2.0)
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)
 
# Track position and momentum
obs = (x = wf -> expectation_value(wf, x -> x), 
       p = wf -> expectation_value(wf, build_momentum_operator(basis)))

t_total = 10.0
dt = 0.05
history, obs_data = propagate_tdse(psi0, ham, t_total, dt; observables=obs)

# Visualization
fig = Figure(size=(800, 600))
ax1 = Axis(fig[1, 1], title="Observables over time", xlabel="Time Step")
lines!(ax1, float.(obs_data.x), label="<x>")
lines!(ax1, real.(obs_data.p), label="<p>")
axislegend(ax1)

ax2 = Axis(fig[2, 1], title="Wavefunction components (Final state)", xlabel="x")
wf_final = history[end]
# Normalize potential for plotting
V_vals = [ham.potential(xi) for xi in basis.x]
V_plot = real.(V_vals) ./ maximum(abs.(V_vals)) .* maximum(abs2.(wf_final.psi))
lines!(ax2, basis.x, V_plot, color=:black, linestyle=:dash, label="V(x) (scaled)")
lines!(ax2, basis.x, real.(wf_final.psi), label="Re(ψ)")
lines!(ax2, basis.x, real.(imag.(wf_final.psi)), label="Im(ψ)")
lines!(ax2, basis.x, real.(abs2.(wf_final.psi)), label="|ψ|²")
axislegend(ax2)

save("examples/results/02_observables.png", fig)
println("Saved results to examples/results/02_observables.png")
