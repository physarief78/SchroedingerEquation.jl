# Basis and Grids

`SchroedingerEquation.jl` uses a real-space discretization for solving the Schrödinger equation. The primary structure for defining the simulation domain is the `RealSpaceGrid1D`.

## RealSpaceGrid1D

The `RealSpaceGrid1D` represents a uniform discretization of a 1D spatial interval $[x_{min}, x_{max}]$ with $N$ points.

### Construction

You can construct a grid using physical units from `Unitful.jl`:

```julia
using SchroedingerEquation
using Unitful

# A grid from -10nm to 10nm with 1000 points
basis = RealSpaceGrid1D(-10.0u"nm", 10.0u"nm", 1000)
```

The constructor automatically converts the input units to internal atomic units (Bohr) for consistent numerical processing.

### Properties

- `basis.x`: A `Vector{Float64}` containing the coordinates of the grid points.
- `basis.dx`: The grid spacing $\Delta x$.
- `basis.N`: The number of points.

## Working with Wavefunctions

Wavefunctions in this package are typically represented by the `Wavefunction1D` struct, which pairs a complex vector of values with its underlying basis.

```julia
# Example: Creating a Gaussian wavepacket on a basis
x = basis.x
psi_vals = exp.(-(x .- 0.0).^2 ./ (2 * 1.0^2))
psi = Wavefunction1D(psi_vals ./ norm(psi_vals), basis)
```

## Momentum Space

For solvers like the **Split-Step Fourier Method (SSFM)**, the package automatically handles the transformation between real space and momentum space ($k$-space).

The momentum grid is defined such that it is compatible with the Fast Fourier Transform (FFT) of the real-space grid.
