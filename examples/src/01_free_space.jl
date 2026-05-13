using SchroedingerEquation
using LinearAlgebra
using CairoMakie

# ==========================================
# 01. Setup the Spatial Grid
# ==========================================
basis = RealSpaceGrid1D(-50.0, 50.0, 1000)

# ==========================================
# 02. Build the Hamiltonian
# ==========================================
# Free particle (V = 0)
potential = CustomPotential(x -> 0.0)
ham = build_hamiltonian(basis, potential)

# ==========================================
# 03. Create Initial Wavepacket
# ==========================================
x0 = -20.0
sigma = 2.0
p0 = 2.0
psi0_vals = [exp(-(x - x0)^2 / (2 * sigma^2)) * exp(im * p0 * x) for x in basis.x]
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)

# ==========================================
# 04. Time Evolution
# ==========================================
t_total = 20.0
dt = 0.02
history = propagate_tdse(psi0, ham, t_total, dt; save_stride=20)

# ==========================================
# 05. Visualize
# ==========================================
fig = Figure(size=(800, 400))
ax = Axis(fig[1,1], title="Free Particle Dispersion", xlabel="x", ylabel="|ψ|²")
lines!(ax, basis.x, abs2.(history[1].psi), label="t=0")
lines!(ax, basis.x, abs2.(history[end].psi), label="t=end")
axislegend(ax)
save("examples/results/01_free_space.png", fig)
println("Saved results to examples/results/01_free_space.png")
