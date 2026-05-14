# SchroedingerEquation.jl

Welcome to the documentation for **SchroedingerEquation.jl**, a high-performance, rigorous, and unit-aware solver for the 1D Schrödinger equation. 

This package was born out of the need for a simulation tool that bridges the gap between educational toy models and hardcore research applications. It is designed for researchers, computational physicists, and students who need robust simulation tools for quantum mechanics, ranging from calculating simple bound states to simulating ultrafast time-dependent laser excitations.

By leveraging the Julia language's speed and ecosystem, `SchroedingerEquation.jl` allows you to define complex potentials and wavepackets and instantly solve for their eigenstates or propagate them through time without sacrificing performance.

## Core Features

- **Strict Unit Support**: Native integration with `Unitful.jl` and `UnitfulAtomic.jl` ensures you never drop a factor of $\hbar$, miscalculate a mass, or get your electron-volts mixed up with your Joules.
- **Fast Solvers**:
    - **TISE (Time-Independent)**: Extremely fast $O(N)$ Eigen-solvers using `SymTridiagonal` matrices and `KrylovKit`. Calculate hundreds of bound states in milliseconds.
    - **TDSE (Time-Dependent)**: Stable Unitary Crank-Nicolson algorithms and ultra-fast Split-Step Fourier Methods (SSFM) for seamlessly propagating wavepackets through time.
- **Advanced Boundaries**: Full support for Periodic boundaries, Hard Walls, and Manolopoulos Complex Absorbing Potentials (CAP) for accurately simulating open quantum systems (like electrons escaping a detector).
- **Visualization**: Beautiful interactive live feedback and professional GIF/MP4 animations seamlessly integrated via `Makie.jl` and `CairoMakie`.

## Getting Started

Let's look at a very simple example: solving for the ground state energy of an electron in a harmonic trap. Notice how we use physical units explicitly!

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

Because the package is currently pending registration in the Julia General Registry, you must install it directly from the GitHub repository URL using the Julia package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/physarief78/SchroedingerEquation.jl")
```
