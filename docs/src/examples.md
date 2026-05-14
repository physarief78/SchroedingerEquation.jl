# Tutorials & Examples

Explore common physics problems solved with `SchroedingerEquation.jl`.

## 1. Free Space Dispersion

This example demonstrates how to set up the spatial grid, initialize a Gaussian wavepacket, and propagate it in free space (V = 0) using the Time-Dependent Schrödinger Equation (TDSE). 

```julia
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
```

**Output:**

![Free Particle Dispersion](assets/01_free_space.png)

As expected, the wavepacket travels to the right and its width increases (disperses) over time!

---

## 2. Quantum Tunneling

In this example, we fire a wavepacket at a rectangular potential barrier. The energy of the wavepacket is lower than the barrier height, so classically it should bounce off. However, quantum mechanics allows it to tunnel through.

```julia
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
animate_evolution(history, barrier; filename="03_tunneling.mp4", fps=30, title="Tunneling Effect")
```

**Output:**

![Quantum Tunneling](assets/03_tunneling.mp4)

Notice how a small fraction of the wavepacket successfully emerges on the right side of the barrier!

---

For more advanced examples, such as multi-threading, python integration, or complex absorbing potentials, check out the `examples/` directory in the repository!
