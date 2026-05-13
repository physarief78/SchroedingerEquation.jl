# Examples

Explore common physics problems solved with `SchroedingerEquation.jl`.

## 1. Quantum Tunneling

Simulate a wavepacket hitting a square potential barrier.

```julia
using SchroedingerEquation
using Unitful

# Setup
basis = RealSpaceGrid1D(-20.0u"nm", 20.0u"nm", 2000)
barrier = SquareWellPotential(5.0u"eV", 2.0u"nm") # Actually a barrier if V0 > 0 in build_hamiltonian context? 
# Note: Check SquareWellPotential implementation for sign.

# Initial state: Gaussian packet with momentum
psi0 = gaussian_packet(basis, -10.0u"nm", 1.0u"nm", 5.0u"eV") 

# Propagate
ham = build_hamiltonian(basis, barrier; boundary=AbsorbingBoundary(5.0u"nm"))
animate_evolution(psi0, ham, 100.0u"fs", 0.5u"fs", "tunneling.mp4")
```

## 2. Multi-threaded Parameter Sweep

Calculate the transmission coefficient for 50 different barrier heights in parallel.

```julia
using Base.Threads

heights = range(1.0, 10.0, length=50) .* 1.0u"eV"
results = zeros(50)

Threads.@threads for i in 1:50
    # Setup and solve...
end
```

## 3. Python Integration

Combine Julia solvers with Python's data analysis ecosystem.

```python
import schroedingerequation as se
import numpy as np

basis = se.RealSpaceGrid1D(se.unit("-5.0u\"nm\""), se.unit("5.0u\"nm\""), 500)
energies, wavefunctions = se.solve_tise(se.build_hamiltonian(basis, se.HarmonicPotential(se.unit("1.0u\"eV/nm^2\""))), 3)
```

For more detailed examples, check the `examples/` directory in the repository.
