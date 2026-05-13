# Python Support

`SchroedingerEquation.jl` provides a high-performance Python interface using [`juliacall`](https://github.com/JuliaPy/PythonCall.jl). This allows you to leverage the speed of Julia's Schrödinger solvers directly from within your Python research environment.

## Installation

To use the Python interface, you can install the package directly from the repository:

```bash
pip install "git+https://github.com/physarief78/SchroedingerEquation.jl.git#egg=schroedingerequation&subdirectory=python"
```

Alternatively, if you have the repository cloned locally:

```bash
cd python
pip install .
```

## Quick Start

The Python API is designed to feel native to Python users, with automatic conversion to and from NumPy arrays.

```python
import schroedingerequation as se
import numpy as np
import matplotlib.pyplot as plt

# 1. Setup a grid with physical units
basis = se.RealSpaceGrid1D(se.unit("-2.0u\"nm\""), se.unit("2.0u\"nm\""), 500)

# 2. Define a potential (e.g., Harmonic Oscillator)
pot = se.HarmonicPotential(se.unit("500.0u\"eV/nm^2\""))

# 3. Solve the Time-Independent Schrödinger Equation (TISE)
ham = se.build_hamiltonian(basis, pot)
energies, psi = se.solve_tise(ham, 3)

# 4. Access results as NumPy arrays
x = np.array(basis.x)
plt.plot(x, psi[:, 0]) # Plot the ground state
plt.show()
```

## Features

### Unit Handling
Physical units are handled via the `se.unit()` helper, which supports both standard physical units (via `Unitful.jl`) and atomic units (via `UnitfulAtomic.jl`).

### NumPy Integration
- **Grid coordinates**: `np.array(basis.x)`
- **Eigenenergies**: Returned as a standard 1D NumPy array.
- **Wavefunctions**: Returned as a 2D NumPy array with shape `(n_points, n_states)`.

## Demo Example
A complete, professional example with high-resolution plotting is available in the repository at `examples/python_demo/src/main.py`.
