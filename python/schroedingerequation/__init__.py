"""
Python wrapper for SchroedingerEquation.jl
"""

import os
from juliacall import Main as jl
import numpy as np

# Find the project root (where Project.toml is)
_base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_project_toml = os.path.join(_base_dir, "Project.toml")

def _init_julia():
    try:
        # Load Pkg and register the local project if we are in the source tree
        jl.eval("import Pkg")
        if os.path.exists(_project_toml):
            path = _base_dir.replace("\\", "/")
            jl.eval(f'Pkg.develop(path="{path}")')
        
        # Load the package into Julia's Main namespace
        jl.eval("using SchroedingerEquation")
        return jl.SchroedingerEquation
    except Exception as e:
        # Fallback: try to add from GitHub if local load fails
        try:
            print(f"SchroedingerEquation.jl local load failed: {e}. Attempting fallback...")
            jl.eval("using SchroedingerEquation")
            return jl.SchroedingerEquation
        except:
            print("Attempting to install SchroedingerEquation.jl from GitHub...")
            jl.eval('Pkg.add(url="https://github.com/physarief78/SchroedingerEquation.jl")')
            jl.eval("using SchroedingerEquation")
            return jl.SchroedingerEquation

se = _init_julia()

# Helper to convert numpy arrays to Julia arrays if needed
def to_jl(x):
    return x

# Export core types and functions from the module
RealSpaceGrid1D = se.RealSpaceGrid1D
Wavefunction1D = se.Wavefunction1D
CustomPotential = se.CustomPotential
CustomTimeDependentPotential = se.CustomTimeDependentPotential
Hamiltonian1D = se.Hamiltonian1D

HarmonicPotential = se.HarmonicPotential
SquareWellPotential = se.SquareWellPotential
OscillatingField = se.OscillatingField
CombinedPotential = se.CombinedPotential
ManolopoulosCAP = se.ManolopoulosCAP

AbsorbingBoundary = se.AbsorbingBoundary
PeriodicBoundary = se.PeriodicBoundary
HardWallBoundary = se.HardWallBoundary

# Functions
probability_density = se.probability_density
expectation_value = se.expectation_value
strip_units = se.strip_units
apply_units = se.apply_units
build_momentum_operator = se.build_momentum_operator
propagate_tdse = se.propagate_tdse
propagate_ssfm = se.propagate_ssfm
animate_evolution = se.animate_evolution

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

# Smoke test
assert hasattr(se, "RealSpaceGrid1D"), "Failed to load core types from SchroedingerEquation"

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
    "probability_density",
    "expectation_value",
    "strip_units",
    "apply_units",
    "build_momentum_operator",
    "propagate_tdse",
    "propagate_ssfm",
    "animate_evolution",
    "build_hamiltonian",
    "solve_tise",
    "simulate_tdse",
    "unit",
]
