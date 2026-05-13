"""
Main entry point for Python demo: Solving TISE for Harmonic Oscillator
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add project root to sys.path to allow importing schroedingerequation without installation
# Current file is in examples/python_demo/src/
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(project_root, "python"))

import schroedingerequation as se

def main():
    print("Initializing SchroedingerEquation from Python...")
    results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))
    os.makedirs(results_dir, exist_ok=True)
    
    # 1. Balanced resolution grid (2nm wide, 800 points)
    basis = se.RealSpaceGrid1D(se.unit("-1.0u\"nm\""), se.unit("1.0u\"nm\""), 800)
    
    # 2. Strong Harmonic Potential (800 eV/nm^2)
    pot = se.HarmonicPotential(se.unit("800.0u\"eV/nm^2\""))
    
    # 3. Build Hamiltonian
    ham = se.build_hamiltonian(basis, pot)
    
    # 4. Solve TISE for first 5 states
    print("Solving TISE...")
    energies, psi = se.solve_tise(ham, 5)
    
    # 5. Visualize with Matplotlib
    x = np.array(basis.x)
    
    plt.figure(figsize=(10, 6))
    for i in range(3):
        # Extract the i-th wavefunction
        wf_raw = np.real(psi[:, i])
        
        # Flip the phase if the peak is negative, so it matches textbook conventions
        peak_idx = np.argmax(np.abs(wf_raw))
        if wf_raw[peak_idx] < 0:
            wf_raw = -wf_raw
            
        # Shift wavefunctions by energy for better visualization
        plt.plot(x, wf_raw + energies[i], label=f"n={i} (E={energies[i]:.4f})")
    
    plt.title("Eigenstates of Harmonic Oscillator (Python Interface)")
    plt.xlabel("Position (Bohr)")
    plt.ylabel("Energy (Hartree) / Amplitude")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Save result
    save_path = os.path.join(results_dir, "eigenstates.png")
    plt.savefig(save_path, dpi=300)
    print(f"Result saved to {save_path}")
    
    # Show plot if interactive
    plt.show()

if __name__ == "__main__":
    main()
