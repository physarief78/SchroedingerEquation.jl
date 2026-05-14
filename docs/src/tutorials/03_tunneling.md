# 3. Quantum Tunneling

In this example, we tackle a phenomenon that is impossible in classical mechanics: Quantum Tunneling. We will fire a wavepacket at a rectangular potential barrier. The kinetic energy of the wavepacket is strictly lower than the barrier's potential height. Classically, it should completely bounce off. Let's see what happens quantum mechanically.

### Setup and Initialization
We use the built-in `SquareWellPotential` to construct a barrier. Notice that by giving it a negative depth (e.g., `-15.0`), it acts as a positive barrier in our Hamiltonian conventions.

```julia
using SchroedingerEquation
using LinearAlgebra
using CairoMakie

# Setup
basis = RealSpaceGrid1D(-50.0, 50.0, 1000)

# A barrier in the center (positive potential height 15)
barrier = SquareWellPotential(2.0, -15.0; x0=0.0) 
ham = build_hamiltonian(basis, barrier)
```

### The Tunneling Regime
We carefully set the momentum `p0 = 5.0`. In our natural units ($m=1, \hbar=1$), the kinetic energy is $E = p^2 / 2m = 12.5$. Since $12.5 < 15.0$, we are firmly in the tunneling regime.

```julia
# Initial state: A moving Gaussian wavepacket starting on the left
x0 = -25.0
sigma = 2.0
p0 = 5.0 
psi0_vals = [exp(-(x - x0)^2 / (2 * sigma^2)) * exp(im * p0 * x) for x in basis.x]
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)

# TDSE: Tunneling animation
t_total = 12.0
dt = 0.02
history = propagate_tdse(psi0, ham, t_total, dt; save_stride=5)
```

### Rendering the Animation
Instead of plotting static lines, we use the incredibly powerful `animate_evolution` tool to generate an MP4 video of the entire process!

```julia
# Visualization
animate_evolution(history, barrier; filename="examples/results/03_tunneling.mp4", fps=30, title="Tunneling Effect")
println("Saved animation to examples/results/03_tunneling.mp4")
```

**Console Output:**
```text
Saved animation to examples/results/03_tunneling.mp4
```

**Resulting Animation:**

![Quantum Tunneling](../assets/03_tunneling.mp4)

Notice how the majority of the wave is reflected, but a small fraction of the wavepacket successfully emerges and continues propagating on the right side of the barrier!
