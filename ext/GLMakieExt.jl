module GLMakieExt

using SchroedingerEquation
using GLMakie
using Makie
using LinearAlgebra
using SparseArrays

function SchroedingerEquation.live_evolution(psi0::Wavefunction1D, ham::Hamiltonian1D, t_total::Float64, dt::Float64; stride=100, hbar=1.0, show_real=true, show_imag=false, show_abs2=true)
    # 1. Setup the GLMakie Figure
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title="Live Quantum Dynamics", xlabel="x", ylabel="ψ(x), V(x)")
    
    # 2. Create Observables for the plot
    basis = ham.basis
    x = basis.x
    wf_obs = Observable(psi0)
    
    # Plot potential in the background
    V_vals = [ham.potential(xi) for xi in x]
    v_max = maximum(abs.(V_vals))
    if v_max > 0
        lines!(ax, x, real.(V_vals) .* (1.0 / v_max) .* 0.5, color=:black, linestyle=:dash, label="V(x)")
    end
    
    # Plot components
    if show_real
        lines!(ax, x, lift(wf -> real.(wf.psi), wf_obs), color=:blue, label="Re(ψ)")
    end
    if show_imag
        lines!(ax, x, lift(wf -> imag.(wf.psi), wf_obs), color=:green, label="Im(ψ)")
    end
    if show_abs2
        lines!(ax, x, lift(wf -> abs2.(wf.psi), wf_obs), color=:red, label="|ψ|²")
    end
    axislegend(ax)
    
    display(fig) # Opens the interactive window
    
    # 3. Fast TDSE Setup (Extract Tridiagonal matrices)
    N = basis.N
    V_matrix = spdiagm(0 => V_vals)
    T_mat = ham.H_matrix - V_matrix
    T_tri = Tridiagonal(T_mat)
    
    A_d = zeros(ComplexF64, N)
    B_d = zeros(ComplexF64, N)
    gamma = im * dt / (2 * hbar)
    A_dl = collect(gamma .* T_tri.dl)
    A_du = collect(gamma .* T_tri.du)
    B_dl = conj.(A_dl)
    B_du = conj.(A_du)
    B_psi = zeros(ComplexF64, N)
    
    psi = copy(psi0.psi)
    n_steps = Int(round(t_total / dt))
    
    # 4. The Live Evolution Loop
    for step in 1:n_steps
        # Check if the user closed the window
        if !isopen(fig.scene)
            println("Live view closed. Stopping simulation.")
            break
        end

        # Crank-Nicolson Math
        @inbounds @simd for j in 1:N
            H_diag = T_tri.d[j] + ham.potential(x[j])
            A_d[j] = 1.0 + gamma * H_diag
            B_d[j] = 1.0 - gamma * H_diag
        end

        A = Tridiagonal(A_dl, A_d, A_du)
        B = Tridiagonal(conj.(A_dl), B_d, conj.(A_du))

        mul!(B_psi, B, psi)
        psi .= A \ B_psi

        # PUSH to the live plot every `stride` steps
        if step % stride == 0
            wf_obs[] = Wavefunction1D(psi, basis)
            sleep(0.001)
        end
    end
end

end # module
