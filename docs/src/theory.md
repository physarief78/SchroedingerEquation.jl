# Theory & Numerical Methods

This document details the underlying physical theory and the numerical methods implemented in **SchroedingerEquation.jl**.

## 1. Physics Theory

### Time-Independent Schrödinger Equation (TISE)
The stationary states of a quantum system are governed by the Time-Independent Schrödinger Equation:

$$\hat{H} \psi_n(x) = E_n \psi_n(x)$$

where $\hat{H}$ is the Hamiltonian operator:

$$\hat{H} = -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} + V(x)$$

$E_n$ is the eigenenergy and $\psi_n(x)$ is the eigenfunction corresponding to the $n$-th state.

### Time-Dependent Schrödinger Equation (TDSE)
The time evolution of a quantum state $\psi(x,t)$ is described by the Time-Dependent Schrödinger Equation:

$$i\hbar \frac{\partial}{\partial t} \psi(x,t) = \hat{H}(t) \psi(x,t)$$

If the Hamiltonian is time-independent, the formal solution is:

$$\psi(x,t) = \exp\left(-\frac{i \hat{H} t}{\hbar}\right) \psi(x,0)$$

## 2. Unit System and Atomic Units

**SchroedingerEquation.jl** natively integrates with `Unitful.jl` and `UnitfulAtomic.jl` to handle physical quantities. All internal computations are strictly performed in **Atomic Units (a.u.)** where $\hbar = 1$, $m_e = 1$, and $e = 1$.

When you provide a parameter with physical units (e.g., `1.0u"eV/nm^2"` or `-10.0u"nm"`), the package automatically uses `UnitfulAtomic.auconvert` to strip the units and scale the values down to their dimensionless atomic representations. After the computation, the user can easily convert the dimensionless results (like energy in Hartrees) back into physical units (like $eV$) if needed.

## 3. Numerical Methods

### Discretization and the TISE
To numerically solve the TISE, space is discretized into $N$ grid points with a uniform spacing $\Delta x$. The kinetic energy operator is approximated using a standard three-point central finite difference scheme:

$$\frac{\partial^2 \psi}{\partial x^2} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\Delta x^2}$$

This maps the continuous operator $\hat{H}$ to a discrete matrix $H$. Because the kinetic energy only couples nearest neighbors and the potential energy is purely diagonal, $H$ takes the form of a **symmetric tridiagonal matrix**. In Julia, we construct this efficiently using `LinearAlgebra.SymTridiagonal`.

The lowest $k$ eigenstates are then found with extreme efficiency in $O(N)$ operations using the `KrylovKit.eigsolve` routine, leveraging the sparse tridiagonal structure.

### Split-Step Fourier Method (SSFM) for TDSE
For simulating the TDSE, especially in large domains or long times, the **Split-Step Fourier Method (SSFM)** is implemented. This method takes advantage of the fact that the kinetic energy operator $\hat{T}$ is diagonal in momentum space, while the potential energy operator $\hat{V}$ is diagonal in position space.

Using the symmetric Strang splitting, the time evolution operator is approximated as:

$$\exp\left(-i \hat{H} \Delta t\right) \approx \exp\left(-i \frac{\hat{V} \Delta t}{2}\right) \exp\left(-i \hat{T} \Delta t\right) \exp\left(-i \frac{\hat{V} \Delta t}{2}\right) + \mathcal{O}(\Delta t^3)$$

The algorithm proceeds as follows for each time step $\Delta t$:
1. Apply half a step of the potential operator in position space.
2. Fourier transform the wavefunction $\psi \xrightarrow{\mathcal{F}} \tilde{\psi}$ using FFT.
3. Apply the full kinetic energy operator in momentum space (multiplying by $\exp(-i \frac{k^2}{2m} \Delta t)$).
4. Inverse Fourier transform back to position space $\tilde{\psi} \xrightarrow{\mathcal{F}^{-1}} \psi$.
5. Apply the second half step of the potential operator.

This spectral method runs in $O(N \log N)$ time, boasts exceptional stability, and inherently imposes periodic boundary conditions.
