import os
import sys

# Add python source to path so we can import without installing
sys.path.insert(0, os.path.abspath("python"))

import numpy as np
import schroedingerequation as se
import matplotlib.pyplot as plt

# 1. Setup grid and potential with physical units
basis = se.RealSpaceGrid1D(se.unit("-5.0u\"nm\""), se.unit("5.0u\"nm\""), 500)
pot = se.HarmonicPotential(se.unit("2.0u\"eV/nm^2\""))

# 2. Build Hamiltonian and solve
ham = se.build_hamiltonian(basis, pot)
energies, wavefunctions = se.solve_tise(ham, 3)

print("Ground state energy:", energies[0])

# 3. Use results in NumPy/Matplotlib
x = np.array(basis.x)
plt.plot(x, np.abs(wavefunctions[:, 0])**2)
plt.savefig("python_output.png")
print("Saved python_output.png")
