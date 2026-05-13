# Potentials and Hamiltonians

The physics of your simulation is defined by the potential energy operator $V(x)$ and the resulting Hamiltonian $H = T + V$.

## Potentials

`SchroedingerEquation.jl` provides several built-in potential types, all of which support physical units.

### Built-in Potentials

- **`HarmonicPotential(k)`**: $V(x) = \frac{1}{2} k x^2$.
- **`SquareWellPotential(V0, width)`**: A finite square well of depth $V0$ and given width.
- **`CustomPotential(f)`**: Define any arbitrary potential using a Julia function `f(x)`.

### Example: Custom Potential

```julia
using SchroedingerEquation
using Unitful

# Define a double well potential
V_func(x) = (x^2 - 5.0^2)^2 * 0.01
pot = CustomPotential(V_func)
```

## Building the Hamiltonian

The `build_hamiltonian` function assembles the kinetic and potential energy operators into a single object.

```julia
ham = build_hamiltonian(basis, pot)
```

### Numerical Implementation

For the Time-Independent Schrödinger Equation (TISE), the Hamiltonian is represented as a **Sparse Matrix** (specifically a `SymTridiagonal` or `SparseMatrixCSC`). This allows for $O(N)$ memory scaling and very fast eigen-solving.

## Boundary Conditions

You can specify boundary conditions when building the Hamiltonian or setting up solvers:

- **`HardWallBoundary()`**: The wavefunction vanishes at the grid edges (Dirichlet).
- **`PeriodicBoundary()`**: The system is treated as periodic (Born-von Karman).
- **`AbsorbingBoundary(L)`**: Adds a Complex Absorbing Potential (CAP) of length $L$ to the edges to simulate open boundaries (scattering).

```julia
# Hamiltonian with Absorbing Boundaries
ham = build_hamiltonian(basis, pot; boundary=AbsorbingBoundary(2.0u"nm"))
```
