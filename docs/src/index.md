# SchroedingerEquation.jl

Welcome to the documentation for **SchroedingerEquation.jl**, a high-performance, rigorous, and unit-aware solver for the 1D Schrödinger equation. 

This package is designed for researchers and students who need robust simulation tools for quantum mechanics, from simple eigenstates to ultrafast time-dependent dynamics.

## Features

- **Unit Support**: Native integration with `Unitful.jl` and `UnitfulAtomic.jl` ensures you never drop a factor of $\hbar$.
- **Fast Solvers**:
    - **TISE (Time-Independent)**: Extremely fast $O(N)$ Eigen-solvers using `SymTridiagonal` matrices and `KrylovKit`.
    - **TDSE (Time-Dependent)**: Unitary Crank-Nicolson and ultra-fast Split-Step Fourier Method (SSFM) for propagating wavepackets.
- **Boundaries**: Support for Periodic, Hard Wall, and Manolopoulos Complex Absorbing Potentials (CAP) for open quantum systems.
- **Visualization**: Interactive live feedback and professional mp4 animations seamlessly integrated via `Makie.jl`.

## Getting Started

Let's solve for the ground state energy of an electron in a harmonic trap:

```julia
using SchroedingerEquation
using Unitful
using UnitfulAtomic

# 1. Define a 20nm spatial grid with 1000 points
basis = RealSpaceGrid1D(-10.0u"nm", 10.0u"nm", 1000)

# 2. Define a harmonic potential trap (k = 1.0 eV/nm²)
pot = HarmonicPotential(1.0u"eV/nm^2")

# 3. Build the Hamiltonian
ham = build_hamiltonian(basis, pot)

# 4. Solve the Time-Independent Schrödinger Equation for the lowest 3 states
energies, states = solve_tise(ham, 3)

# 5. Extract and print the ground state energy in eV
println("Ground state energy: ", apply_units(energies[1], u"eV"))
```

**Output:**
```text
Ground state energy: 4.532029953521174e-5 eV
```

## Installation

You can install `SchroedingerEquation.jl` directly via the Julia package manager:

```julia
using Pkg
Pkg.add("SchroedingerEquation")
```
