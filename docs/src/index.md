# SchroedingerEquation.jl

Welcome to the documentation for **SchroedingerEquation.jl**, a high-performance, rigorous, and unit-aware solver for the 1D Schrödinger equation.

## Features

- **Unit Support**: Native integration with `Unitful.jl` and `UnitfulAtomic.jl`.
- **Fast Solvers**:
    - **TISE**: $O(N)$ Eigen-solvers using `SymTridiagonal` matrices and `KrylovKit`.
    - **TDSE**: Unitary Crank-Nicolson and ultra-fast Split-Step Fourier Method (SSFM).
- **Boundaries**: Periodic, Hard Wall, and Manolopoulos Complex Absorbing Potentials (CAP).
- **Visualization**: Interactive live feedback and professional animations using `Makie.jl`.

## Quick Start

```julia
using SchroedingerEquation
using Unitful

# 1. Define Basis (Unit-aware!)
basis = RealSpaceGrid1D(-10.0u"nm", 10.0u"nm", 1000)

# 2. Define Potential
pot = HarmonicPotential(1.0u"eV/nm^2")

# 3. Build Hamiltonian
ham = build_hamiltonian(basis, pot)

# 4. Solve TISE
energies, states = solve_tise(ham, 3)

println("Ground state energy: ", energies[1])
```

## Installation

```julia
using Pkg
Pkg.add("SchroedingerEquation")
```
