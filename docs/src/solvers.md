# Solvers

`SchroedingerEquation.jl` provides optimized solvers for both the Time-Independent (TISE) and Time-Dependent (TDSE) Schrödinger Equation.

## Time-Independent (TISE)

Solve $H\psi = E\psi$ to find the stationary states (eigenstates) of the system.

```julia
# Find the 5 lowest energy states
energies, wavefunctions = solve_tise(ham, 5)
```

The TISE solver uses `KrylovKit.jl` for high-performance iterative diagonalization, making it suitable for very large grids.

## Time-Dependent (TDSE)

Propagate a wavefunction $\psi(t)$ according to $i\hbar \frac{\partial}{\partial t}\psi = H\psi$.

### 1. Unitary Crank-Nicolson

A stable, 2nd-order implicit method that preserves the norm of the wavefunction. Best for general potentials and non-periodic boundaries.

```julia
psi_final = propagate_tdse(psi_init, ham, t_total, dt)
```

### 2. Split-Step Fourier Method (SSFM)

An ultra-fast spectral method that uses FFTs to solve the kinetic and potential parts separately. It is $O(N \log N)$ and highly accurate for smooth potentials and periodic boundaries.

```julia
psi_final = propagate_ssfm(psi_init, pot, t_total, dt)
```

## Performance Tips

### Workspaces for Parallelism

To avoid redundant allocations in loops (e.g., parameter sweeps), use the **Workspace Pattern**:

```julia
# Pre-allocate workspace
ws = SSFMWorkspace(psi_init, pot, dt)

# Use workspace in a loop
for i in 1:100
    propagate_ssfm(psi_init, potentials[i], t_total, dt; workspace=ws)
end
```

This approach is highly recommended for multi-threaded simulations using `Threads.@threads`.
