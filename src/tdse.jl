# Time-Dependent Schroedinger Equation (TDSE) solvers

"""
    step_crank_nicolson(psi, H_matrix, dt; hbar=1.0)
Perform one time step using the Crank-Nicolson method.
"""
function step_crank_nicolson(psi::Vector{ComplexF64}, H_matrix::AbstractMatrix, dt::Float64; hbar=1.0)
    # Solve A * psi_next = B * psi
    # Using I (UniformScaling) with sparse/dense matrices is more efficient than I(N)
    A = I + (im * dt / (2 * hbar)) * H_matrix
    B = I - (im * dt / (2 * hbar)) * H_matrix
    
    psi_next = A \ (B * psi)
    return psi_next
end

"""
    TDSEWorkspace(ham::Hamiltonian1D, dt; hbar=1.0)
Pre-allocates buffers and factorizes the Crank-Nicolson matrix for a static Hamiltonian.
"""
function TDSEWorkspace(ham::Hamiltonian1D, dt; hbar=1.0)
    _hbar = strip_units(hbar)
    _dt = strip_units(dt)
    H = ham.H_matrix
    N = size(H, 1)
    
    gamma = im * _dt / (2 * _hbar)
    A = I + gamma * H
    B = I - gamma * H
    
    return TDSEWorkspace(lu(A), B, zeros(ComplexF64, N))
end

"""
    propagate_tdse(psi0::Wavefunction1D, ham::Hamiltonian1D, t_total, dt; hbar=1.0, save_stride=1, observables=nothing, keep_history=true, workspace=nothing)
Propagate the wavefunction in time for a static Hamiltonian.
Uses pre-factorization of the Crank-Nicolson matrix for speed.
If a `workspace::TDSEWorkspace` is provided, it avoids all internal allocations.
"""
function propagate_tdse(psi0::Wavefunction1D, ham::Hamiltonian1D, t_total, dt; hbar=1.0, save_stride=1, observables=nothing, keep_history=true, workspace=nothing)
    # Handle Unitful inputs
    _hbar = strip_units(hbar)
    _dt = strip_units(dt)
    _t_total = strip_units(t_total)
    
    n_steps = Int(round(_t_total / _dt))
    psi = copy(psi0.psi)
    basis = ham.basis
    N = length(psi)
    
    # Use workspace if provided, otherwise create a temporary one
    ws = if workspace !== nothing
        workspace
    else
        TDSEWorkspace(ham, _dt; hbar=_hbar)
    end
    
    A_fact = ws.A_fact
    B = ws.B
    B_psi = ws.B_psi
    
    history = keep_history ? [Wavefunction1D(copy(psi), basis)] : nothing
    
    # Initialize observables history using NamedTuple and map for type stability
    if observables !== nothing
        obs_history = map(f -> [f(psi0)], observables)
    else
        obs_history = nothing
    end
    
    for i in 1:n_steps
        # Fast in-place multiplication for B * psi
        mul!(B_psi, B, psi)
        psi .= A_fact \ B_psi
        
        if keep_history && (i % save_stride == 0)
            push!(history, Wavefunction1D(copy(psi), basis))
        end
        
        if observables !== nothing
            wf_current = Wavefunction1D(psi, basis)
            for (name, f) in pairs(observables)
                push!(obs_history[name], f(wf_current))
            end
        end
    end
    
    if observables !== nothing
        return keep_history ? (history, obs_history) : obs_history
    else
        return history
    end
end

"""
    propagate_tdse(psi0::Wavefunction1D, basis::RealSpaceGrid1D, T_matrix::AbstractMatrix, potential::AbstractTimeDependentPotential, t_total::Float64, dt::Float64; hbar=1.0, save_stride=1, observables=nothing, keep_history=true)
Propagate the wavefunction in time for a time-dependent potential.
Evaluates the potential at the midpoint of each timestep.
Uses optimized Tridiagonal in-place updates for zero allocation in the loop.
"""
function propagate_tdse(psi0::Wavefunction1D, basis::RealSpaceGrid1D, T_matrix::AbstractMatrix, potential::AbstractTimeDependentPotential, t_total, dt; hbar=1.0, save_stride=1, observables=nothing, keep_history=true)
    # Handle Unitful inputs
    _hbar = strip_units(hbar)
    _dt = strip_units(dt)
    _t_total = strip_units(t_total)
    
    n_steps = Int(round(_t_total / _dt))
    psi = copy(psi0.psi)
    N = basis.N
    
    # Pre-calculate static parts of Crank-Nicolson operators
    T_tri = T_matrix isa Tridiagonal ? T_matrix : Tridiagonal(T_matrix)
    gamma = im * _dt / (2 * _hbar)
    
    # Pre-allocate and pre-calculate components (Zero allocation in loop)
    A_dl = collect(gamma .* T_tri.dl)
    A_du = collect(gamma .* T_tri.du)
    B_dl = conj.(A_dl)
    B_du = conj.(A_du)
    
    A_d = zeros(ComplexF64, N)
    B_d = zeros(ComplexF64, N)
    B_psi = zeros(ComplexF64, N)
    
    history = keep_history ? [Wavefunction1D(copy(psi), basis)] : nothing
    
    # Initialize observables history using NamedTuple and map for type stability
    if observables !== nothing
        obs_history = map(f -> [f(psi0)], observables)
    else
        obs_history = nothing
    end
    
    t = 0.0
    for i in 1:n_steps
        t_mid = t + _dt / 2.0
        
        # 1. Update diagonals in-place (Zero allocation + SIMD Vectorized)
        @inbounds @simd for j in 1:N
            H_diag = T_tri.d[j] + strip_units(potential(basis.x[j], t_mid))
            A_d[j] = 1.0 + gamma * H_diag
            B_d[j] = 1.0 - gamma * H_diag
        end
        
        # 2. Construct lightweight Tridiagonal wrappers (Zero allocation)
        A = Tridiagonal(A_dl, A_d, A_du)
        B = Tridiagonal(B_dl, B_d, B_du)
        
        # 3. Fast Thomas-Algorithm solve
        mul!(B_psi, B, psi)
        psi .= A \ B_psi
        
        if keep_history && (i % save_stride == 0)
            push!(history, Wavefunction1D(copy(psi), basis))
        end
        
        if observables !== nothing
            wf_current = Wavefunction1D(psi, basis)
            for (name, f) in pairs(observables)
                push!(obs_history[name], f(wf_current))
            end
        end
        
        t += _dt
    end
    
    if observables !== nothing
        return keep_history ? (history, obs_history) : obs_history
    else
        return history
    end
end

"""
    SSFMWorkspace(psi0::Wavefunction1D, potential::AbstractPotential, dt; hbar=1.0, m=1.0)
Pre-allocates FFT plans and momentum-space operators for the Split-Step Fourier Method.
"""
function SSFMWorkspace(psi0::Wavefunction1D, potential::AbstractPotential, dt; hbar=1.0, m=1.0)
    _hbar = strip_units(hbar)
    _m = strip_units(m)
    _dt = strip_units(dt)
    
    basis = psi0.basis
    N = basis.N
    dx = basis.dx
    
    # Momentum grid
    k = fftfreq(N, 1.0/dx) .* 2π
    
    # Pre-calculate operators
    V_vals = [strip_units(potential(xi)) for xi in basis.x]
    U_V = exp.(-im * V_vals .* (_dt / (2 * _hbar)))
    
    E_k = (_hbar .* k).^2 ./ (2 * _m)
    U_T = exp.(-im * E_k .* (_dt / _hbar))
    
    # Create plans
    tmp_psi = copy(psi0.psi)
    P_plan = plan_fft(tmp_psi)
    P_inv_plan = plan_ifft(tmp_psi)
    
    return SSFMWorkspace(P_plan, P_inv_plan, U_V, U_T)
end

"""
    propagate_ssfm(psi0::Wavefunction1D, potential::AbstractPotential, t_total, dt; hbar=1.0, m=1.0, save_stride=1, observables=nothing, keep_history=true, workspace=nothing)
Propagate the wavefunction using the Split-Step Fourier Method (SSFM).
Extremely fast for smooth potentials and periodic boundaries. O(N log N).
If a `workspace::SSFMWorkspace` is provided, it avoids plan creation and pre-calculation in every call.
"""
function propagate_ssfm(psi0::Wavefunction1D, potential::AbstractPotential, t_total, dt; hbar=1.0, m=1.0, save_stride=1, observables=nothing, keep_history=true, workspace=nothing)
    # Handle units
    _hbar = strip_units(hbar)
    _m = strip_units(m)
    _dt = strip_units(dt)
    _t_total = strip_units(t_total)
    
    n_steps = Int(round(_t_total / _dt))
    basis = psi0.basis
    psi = copy(psi0.psi)
    
    # Use workspace or create temporary
    ws = if workspace !== nothing
        workspace
    else
        SSFMWorkspace(psi0, potential, _dt; hbar=_hbar, m=_m)
    end
    
    U_V = ws.U_V
    U_T = ws.U_T
    P = ws.P_plan
    P_inv = ws.P_inv_plan
    
    history = keep_history ? [Wavefunction1D(copy(psi), basis)] : nothing
    if observables !== nothing
        obs_history = map(f -> [f(psi0)], observables)
    else
        obs_history = nothing
    end

    for i in 1:n_steps
        # 1. Half-step Potential
        psi .*= U_V
        
        # 2. Full-step Kinetic (in Momentum Space)
        psi_k = P * psi
        psi_k .*= U_T
        psi = P_inv * psi_k
        
        # 3. Half-step Potential
        psi .*= U_V
        
        if keep_history && (i % save_stride == 0)
            push!(history, Wavefunction1D(copy(psi), basis))
        end
        
        if observables !== nothing
            wf_current = Wavefunction1D(psi, basis)
            for (name, f) in pairs(observables)
                push!(obs_history[name], f(wf_current))
            end
        end
    end
    
    if observables !== nothing
        return keep_history ? (history, obs_history) : obs_history
    else
        return history
    end
end

export propagate_tdse, propagate_ssfm
