# Interactive Scattering

One of the standout features of `SchroedingerEquation.jl` is its ability to render **interactive, live-updating** simulations using `GLMakie`.

Because this generates a graphical window on your local machine, we cannot render the output directly in this documentation. However, we encourage you to copy this code and run it in your local Julia REPL!

### Setup
We set up a simple free-space grid with an initial Gaussian wavepacket.

```julia
using SchroedingerEquation
using LinearAlgebra
using GLMakie # GLMakie is required for live_evolution

# Setup
basis = RealSpaceGrid1D(-50.0, 50.0, 1000)
potential = CustomPotential(x -> 0.0)
ham = build_hamiltonian(basis, potential)

# Initial state
x0 = -20.0
sigma = 2.0
p0 = 2.0
psi0_vals = [exp(-(x - x0)^2 / (2 * sigma^2)) * exp(im * p0 * x) for x in basis.x]
psi0 = Wavefunction1D(ComplexF64.(psi0_vals), basis)
normalize!(psi0)
```

### Live Evolution
By calling `live_evolution`, a high-performance OpenGL window will pop up. The simulation will run in a separate thread, continuously updating the wavefunction plot in real-time.

```julia
# Live interactive visualization
println("Starting interactive simulation...")
live_evolution(psi0, ham, 0.02; update_step=5)
```

Try modifying the potential or the starting momentum `p0` to watch how the quantum dynamics change before your eyes!
