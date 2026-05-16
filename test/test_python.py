"""
Tests for Python interface of SchroedingerEquation.jl
"""

import unittest
import numpy as np
import os
import sys

# Add the project root to sys.path so we can import schroedingerequation without installing
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "python")))

import schroedingerequation as se

class TestPythonInterface(unittest.TestCase):
    def test_grid_creation(self):
        # Create a grid
        basis = se.RealSpaceGrid1D(se.unit("-5.0u\"nm\""), se.unit("5.0u\"nm\""), 100)
        self.assertEqual(len(basis.x), 100)
        self.assertIsInstance(np.array(basis.x), np.ndarray)

    def test_tise_solver(self):
        # Harmonic oscillator
        basis = se.RealSpaceGrid1D(se.unit("-10.0u\"nm\""), se.unit("10.0u\"nm\""), 500)
        pot = se.HarmonicPotential(se.unit("1.0u\"eV/nm^2\""))
        ham = se.build_hamiltonian(basis, pot)
        
        # Solve
        energies, psi = se.solve_tise(ham, 3)
        
        # Check results
        self.assertEqual(len(energies), 3)
        self.assertEqual(psi.shape, (500, 3))
        # Ground state energy should be positive
        self.assertGreater(energies[0], 0)
        # First excited state should be higher energy
        self.assertGreater(energies[1], energies[0])

    def test_tise_solver_controls(self):
        basis = se.RealSpaceGrid1D(se.unit("-10.0u\"nm\""), se.unit("10.0u\"nm\""), 500)
        pot = se.HarmonicPotential(se.unit("1.0u\"eV/nm^2\""))
        ham = se.build_hamiltonian(basis, pot)

        energies, psi = se.solve_tise(ham, 3, tol=1e-6, maxiter=300, krylovdim=30)
        self.assertEqual(len(energies), 3)
        self.assertEqual(psi.shape, (500, 3))

        energies, psi, info = se.solve_tise(
            ham,
            3,
            tol=1e-6,
            maxiter=300,
            krylovdim=30,
            return_info=True,
        )
        self.assertEqual(len(energies), 3)
        self.assertEqual(psi.shape, (500, 3))
        self.assertGreaterEqual(int(info.converged), 3)
        self.assertGreaterEqual(len(info.normres), 3)
        self.assertGreater(int(info.numops), 0)

    def test_units(self):
        # Test unit helper
        u = se.unit("1.0u\"eV\"")
        # In atomic units, 1 eV is approx 0.0367
        val = se.strip_units(u)
        self.assertAlmostEqual(float(val), 0.0367493, places=5)

if __name__ == "__main__":
    unittest.main()
