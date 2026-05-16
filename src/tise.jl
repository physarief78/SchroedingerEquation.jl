# Time-Independent Schroedinger Equation (TISE) solvers

"""
    build_hamiltonian(basis, potential, boundary::AbstractBoundaryCondition; hbar=1.0, m=1.0)
Commercial-grade builder that automatically bakes boundary conditions into the static matrix.
"""
function build_hamiltonian(basis::RealSpaceGrid1D, 
                           potential::AbstractPotential, 
                           boundary::AbstractBoundaryCondition=PeriodicBoundary(); 
                           hbar=1.0, m=1.0)
    N = basis.N
    dx = basis.dx
    x = basis.x

    # Handle Unitful inputs
    _hbar = strip_units(hbar)
    _m = strip_units(m)

    # Kinetic Energy Matrix
    coeff = -_hbar^2 / (2 * _m * dx^2)
    
    # 1. Build basic Hamiltonian (Tridiagonal)
    # We use SymTridiagonal if possible for O(N) performance in eigen-solvers
    main_diag = fill(-2 * coeff, N)
    off_diag = fill(coeff, N - 1)
    
    # Add potential to main diagonal
    for i in 1:N
        main_diag[i] += strip_units(potential(x[i]))
    end

    # 2. Apply Boundary Conditions and CAP
    if boundary isa PeriodicBoundary
        # Periodic boundaries break Tridiagonal structure, must use Sparse
        H_matrix = spdiagm(-1 => off_diag, 0 => main_diag, 1 => off_diag)
        H_matrix[1, N] = H_matrix[N, 1] = coeff
    elseif boundary isa AbsorbingBoundary
        # CAP adds complex imaginary parts, must use Sparse
        H_matrix = spdiagm(-1 => off_diag, 0 => main_diag, 1 => off_diag)
        cap = ManolopoulosCAP(basis, boundary; m=_m, hbar=_hbar)
        V_cap = spdiagm(0 => cap.grid_profile)
        H_matrix += V_cap
    else
        # Standard HardWall boundary preserves SymTridiagonal!
        # This is the high-performance path for TISE
        H_matrix = SymTridiagonal(main_diag, off_diag)
    end
    
    return Hamiltonian1D(basis, potential, H_matrix)
end

"""
    build_momentum_operator(basis::RealSpaceGrid1D; hbar=1.0)
Constructs the momentum operator matrix using central finite difference.
"""
function build_momentum_operator(basis::RealSpaceGrid1D; hbar=1.0)
    N = basis.N
    dx = basis.dx

    # First derivative central difference: f'(x) ≈ (f(x+dx) - f(x-dx)) / 2dx
    coeff = -im * hbar / (2 * dx)
    off_diag = fill(1.0, N - 1)
    
    # +1 on upper diagonal, -1 on lower diagonal
    P_matrix = coeff * spdiagm(-1 => -off_diag, 1 => off_diag)
    
    return P_matrix
end

function _validate_tise_controls(n_states::Int, tol::Real, maxiter::Int, krylovdim::Int)
    n_states > 0 || throw(ArgumentError("n_states must be positive, got $n_states"))
    tol > 0 || throw(ArgumentError("tol must be positive, got $tol"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive, got $maxiter"))
    krylovdim >= n_states || throw(ArgumentError("krylovdim must be at least n_states ($n_states), got $krylovdim"))
    return nothing
end

"""
    solve_tise(hamiltonian::Hamiltonian1D, n_states::Int; tol=KrylovDefaults.tol[], maxiter=KrylovDefaults.maxiter[], krylovdim=KrylovDefaults.krylovdim[], verbosity=KrylovDefaults.verbosity[], return_info=false)
Solves the TISE for the lowest `n_states` energy levels using KrylovKit.
Returns eigenvalues (energies) and eigenvectors (wavefunctions). If
`return_info=true`, also returns KrylovKit's convergence information.

The Krylov solver stops when the residual norm reaches `tol`, rebuilding a
Krylov subspace of size `krylovdim` up to `maxiter` times.
"""
function solve_tise(hamiltonian::Hamiltonian1D, n_states::Int;
                    tol::Real=KrylovDefaults.tol[],
                    maxiter::Int=KrylovDefaults.maxiter[],
                    krylovdim::Int=KrylovDefaults.krylovdim[],
                    verbosity::Int=KrylovDefaults.verbosity[],
                    return_info::Bool=false)
    _validate_tise_controls(n_states, tol, maxiter, krylovdim)

    H = hamiltonian.H_matrix
    # Use KrylovKit to find the smallest eigenvalues
    vals, vecs, info = eigsolve(H, n_states, :SR; tol=tol, maxiter=maxiter, krylovdim=krylovdim, verbosity=verbosity) # :SR = Smallest Real

    if info.converged < n_states
        @warn "TISE eigensolver did not converge all requested states" requested=n_states converged=info.converged tol=tol normres=info.normres numiter=info.numiter numops=info.numops
    end
    
    energies = real.(vals[1:n_states])
    wavefunctions = [Wavefunction1D(vecs[i], hamiltonian.basis) for i in 1:n_states]
    
    # Ensure wavefunctions are normalized with the dx factor for physical expectation values
    for wf in wavefunctions
        normalize!(wf)
    end
    
    return return_info ? (energies, wavefunctions, info) : (energies, wavefunctions)
end

export build_hamiltonian, build_momentum_operator, solve_tise
