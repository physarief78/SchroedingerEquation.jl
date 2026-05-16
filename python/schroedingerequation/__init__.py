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
    """Load SchroedingerEquation.jl and return the module object."""
    # juliapkg.json handles Pkg.develop automatically via "dev": true.
    # We just need to call 'using' to load the module into Main.
    try:
        jl.seval("using SchroedingerEquation")
        return jl.SchroedingerEquation
    except Exception as e:
        # Fallback: try to install from GitHub
        print(f"SchroedingerEquation.jl failed to load: {e}")
        print("Attempting to install SchroedingerEquation.jl from GitHub...")
        jl.seval("import Pkg")
        jl.seval('Pkg.add(url="https://github.com/physarief78/SchroedingerEquation.jl")')
        jl.seval("using SchroedingerEquation")
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

def solve_tise(
    ham,
    n_states,
    tol=None,
    maxiter=None,
    krylovdim=None,
    verbosity=None,
    return_info=False,
):
    """
    Solves the Time-Independent Schroedinger Equation with optional Krylov controls.
    Returns (energies, wavefunctions) as numpy-friendly objects.
    If return_info=True, also returns KrylovKit convergence information.
    """
    kwargs = {}
    if tol is not None: kwargs["tol"] = tol
    if maxiter is not None: kwargs["maxiter"] = maxiter
    if krylovdim is not None: kwargs["krylovdim"] = krylovdim
    if verbosity is not None: kwargs["verbosity"] = verbosity
    if return_info: kwargs["return_info"] = True

    result = se.solve_tise(ham, n_states, **kwargs)
    if return_info:
        energies, wavefunctions, info = result
    else:
        energies, wavefunctions = result
    
    # Extract raw psi vectors from the list of Wavefunction1D objects
    # and stack them into a 2D NumPy array (points, states)
    psi_list = [np.array(wf.psi) for wf in wavefunctions]
    psi_matrix = np.column_stack(psi_list)
    
    if return_info:
        return np.array(energies), psi_matrix, info
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
    return jl.SchroedingerEquation.eval(jl.Meta.parse(u_str))

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
