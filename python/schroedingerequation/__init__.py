"""
Python wrapper for SchroedingerEquation.jl
"""

import os
from juliacall import Main as jl
import numpy as np

# Initialize Julia environment
try:
    jl.eval(jl.Meta.parse("using SchroedingerEquation"))
    se = jl.eval(jl.Meta.parse("SchroedingerEquation"))
except Exception:
    # If it fails, try to install it from the repository
    jl.eval(jl.Meta.parse("import Pkg"))
    
    # Try to see if we are in a local dev environment first
    project_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if os.path.exists(os.path.join(project_path, "Project.toml")):
        jl.Pkg.activate(project_path)
        jl.Pkg.instantiate()
    else:
        # Otherwise, install from GitHub
        print("SchroedingerEquation.jl not found. Attempting to install from GitHub...")
        jl.Pkg.add(url="https://github.com/physarief78/SchroedingerEquation.jl")
    
    jl.eval(jl.Meta.parse("using SchroedingerEquation"))
    se = jl.eval(jl.Meta.parse("SchroedingerEquation"))

# Helper to convert numpy arrays to Julia arrays if needed
def to_jl(x):
    return x

# Export core types and functions
RealSpaceGrid1D = se.RealSpaceGrid1D
HarmonicPotential = se.HarmonicPotential
SquareWellPotential = se.SquareWellPotential
OscillatingField = se.OscillatingField
CombinedPotential = se.CombinedPotential
ManolopoulosCAP = se.ManolopoulosCAP
AbsorbingBoundary = se.AbsorbingBoundary
PeriodicBoundary = se.PeriodicBoundary
HardWallBoundary = se.HardWallBoundary

def build_hamiltonian(basis, potential, hbar=None, m=None):
    """
    Wrap build_hamiltonian with unit handling if strings are passed.
    """
    kwargs = {}
    if hbar is not None: kwargs['hbar'] = hbar
    if m is not None: kwargs['m'] = m
    return se.build_hamiltonian(basis, potential, **kwargs)

def solve_tise(ham, n_states):
    """
    Solves the Time-Independent Schroedinger Equation.
    Returns (energies, wavefunctions) as numpy-friendly objects.
    """
    energies, wavefunctions = se.solve_tise(ham, n_states)
    
    # Extract raw psi vectors from the list of Wavefunction1D objects
    # and stack them into a 2D NumPy array (points, states)
    psi_list = [np.array(wf.psi) for wf in wavefunctions]
    psi_matrix = np.column_stack(psi_list)
    
    return np.array(energies), psi_matrix

def simulate_tdse(psi0, pot, dt, t_total, method="ssfm", **kwargs):
    """
    Simulates the Time-Dependent Schroedinger Equation.
    """
    return se.simulate_tdse(psi0, pot, dt, t_total, method=method, **kwargs)

# Unit helper
def unit(u_str):
    """Helper to create Unitful units from strings."""
    # Ensure Unitful and UnitfulAtomic are loaded in Main
    jl.eval(jl.Meta.parse("using Unitful, UnitfulAtomic"))
    return jl.eval(jl.Meta.parse(u_str))

__all__ = [
    "RealSpaceGrid1D",
    "HarmonicPotential",
    "SquareWellPotential",
    "OscillatingField",
    "CombinedPotential",
    "ManolopoulosCAP",
    "AbsorbingBoundary",
    "PeriodicBoundary",
    "HardWallBoundary",
    "build_hamiltonian",
    "solve_tise",
    "simulate_tdse",
    "unit",
]
