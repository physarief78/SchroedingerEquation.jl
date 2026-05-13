using Makie
using CairoMakie

# --- Recipes ---

@recipe(PlotWavefunction, wf) do scene
    Attributes(
        color = :blue,
        linewidth = 2
    )
end

function Makie.plot!(p::PlotWavefunction{<:Tuple{Wavefunction1D}})
    wf = p[:wf][]
    x = wf.basis.x
    psi = wf.psi
    
    # Plot real part and magnitude squared
    lines!(p, x, real.(psi), color=p[:color], linewidth=p[:linewidth], label="Re(ψ)")
    lines!(p, x, abs2.(psi), color=:red, linewidth=p[:linewidth], label="|ψ|²")
    
    return p
end

"""
    animate_evolution(history::Vector{Wavefunction1D}, potential=nothing; filename="evolution.mp4", fps=30, title="Wavefunction Evolution", show_real=true, show_imag=false, show_abs2=true)
Generates an animation of the wavefunction evolution from a history array.
If a potential (AbstractPotential, function, or time-dependent potential) is provided, it's plotted in the background.
"""
function animate_evolution(history::Vector{Wavefunction1D}, potential=nothing; filename="evolution.mp4", fps=30, title="Wavefunction Evolution", show_real=true, show_imag=false, show_abs2=true)
    basis = history[1].basis
    x = basis.x
    
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title=title, xlabel="x", ylabel="ψ(x), V(x)")
    
    # Observable for the current wavefunction and time
    wf_obs = Observable(history[1])
    idx_obs = Observable(1)
    
    # Plot potential if provided
    if potential !== nothing
        # We might need to handle time-dependent potentials
        if potential isa AbstractTimeDependentPotential || (potential isa Function && methods(potential)[1].nargs == 3) # (p, x, t)
            v_plot = lift(idx -> begin
                # Note: history doesn't track t, so we use a dummy or first step t
                [potential(xi, 0.0) for xi in x]
            end, idx_obs)
        else
            v_plot = Observable([potential(xi) for xi in x])
        end
        
        # --- NEW ANCHORING LOGIC ---
        # 1. Find the maximum height of the initial wavepacket for visual scaling
        max_psi2 = maximum(abs2.(history[1].psi))
        
        # 2. Scale potential for visualization
        v_vals = v_plot[]
        v_max = maximum(abs.(v_vals))
        
        if v_max > 0
            # Anchor baseline at 0.0, scale peak to 1.2x the height of the wavepacket
            v_scaled = lift(v -> real.(v) .* (max_psi2 / v_max) .* 1.2, v_plot)
            lines!(ax, x, v_scaled, color=:black, linestyle=:dash, label="V(x)", linewidth=2.0)
        end
        # ---------------------------
    end

    # Plot components based on config
    if show_real
        lines!(ax, x, lift(wf -> real.(wf.psi), wf_obs), color=:blue, linewidth=2, label="Re(ψ)")
    end
    if show_imag
        lines!(ax, x, lift(wf -> imag.(wf.psi), wf_obs), color=:green, linewidth=2, label="Im(ψ)")
    end
    if show_abs2
        lines!(ax, x, lift(wf -> abs2.(wf.psi), wf_obs), color=:red, linewidth=2, label="|ψ|²")
    end
    
    axislegend(ax)
    
    # Fix limits
    ylims!(ax, -1.1, 1.1)
    
    record(fig, filename, 1:length(history); framerate=fps) do i
        idx_obs[] = i
        wf_obs[] = history[i]
    end
    
    return filename
end

export animate_evolution
