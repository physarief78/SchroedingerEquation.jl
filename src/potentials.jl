# Common physical potentials

"""
    HarmonicPotential(omega; m=1.0, x0=0.0)
V(x) = 1/2 * m * omega^2 * (x - x0)^2
"""
struct HarmonicPotential <: AbstractPotential
    omega::Float64
    m::Float64
    x0::Float64

    function HarmonicPotential(omega; m=1.0, x0=0.0)
        new(strip_units(omega), strip_units(m), strip_units(x0))
    end
end

(p::HarmonicPotential)(x) = 0.5 * p.m * p.omega^2 * (x - p.x0)^2

"""
    SquareWellPotential(width, depth; x0=0.0)
V(x) = -depth if |x - x0| < width/2, else 0
"""
struct SquareWellPotential <: AbstractPotential
    width::Float64
    depth::Float64
    x0::Float64

    function SquareWellPotential(width, depth; x0=0.0)
        new(strip_units(width), strip_units(depth), strip_units(x0))
    end
end

(p::SquareWellPotential)(x) = abs(x - p.x0) < p.width/2 ? -p.depth : 0.0

"""
    OscillatingField(amplitude, frequency)
Represents a continuous laser field interacting via the dipole approximation: V(x,t) = -x * A * cos(omega * t)
"""
struct OscillatingField <: AbstractTimeDependentPotential
    amplitude::Float64
    omega::Float64

    function OscillatingField(amplitude, frequency)
        new(strip_units(amplitude), strip_units(frequency))
    end
end

(p::OscillatingField)(x, t) = -x * p.amplitude * cos(p.omega * t)

"""
    CombinedPotential(static_pot, dynamic_pot)
Allows adding a static background potential to a time-dependent one.
"""
struct CombinedPotential <: AbstractTimeDependentPotential
    static_pot::AbstractPotential
    dynamic_pot::AbstractTimeDependentPotential
end

(p::CombinedPotential)(x, t) = p.static_pot(x) + p.dynamic_pot(x, t)

"""
    ManolopoulosCAP 
    
Holds the pre-computed grid profile for the absorbing potential.
"""
struct ManolopoulosCAP 
    grid_profile::Vector{ComplexF64}
end

"""
    ManolopoulosCAP(basis::RealSpaceGrid1D, boundary::AbsorbingBoundary; m=1.0, hbar=1.0)
Auto-generates a commercial-grade absorbing boundary. The strength is automatically 
derived from the fundamental constants and grid spacing.
"""
function ManolopoulosCAP(basis::RealSpaceGrid1D, boundary::AbsorbingBoundary; m=1.0, hbar=1.0)
    profile = zeros(ComplexF64, basis.N)
    L = boundary.L
    
    # The pre-factor based on grid physics
    # We include the strength multiplier for fine-tuning if needed (default 1.0)
    coeff = -im * boundary.strength * (hbar^2 / (2.0 * m * L^2))
    
    for i in 1:basis.N
        x = basis.x[i]
        
        # Left boundary
        if x < basis.x_min + L
            z = (basis.x_min + L - x) / L
            profile[i] = coeff * 4.0 * (1.0/(1.0 - z)^2 + 1.0/(1.0 + z)^2 - 2.0)
            
        # Right boundary
        elseif x > basis.x_max - L
            z = (x - (basis.x_max - L)) / L
            profile[i] = coeff * 4.0 * (1.0/(1.0 - z)^2 + 1.0/(1.0 + z)^2 - 2.0)
        end
    end
    
    # Cap the infinity at the absolute edge to prevent solver failures
    # We replace non-finite values (Inf at the boundary) with the maximum finite value
    finite_indices = isfinite.(profile)
    if any(finite_indices)
        max_val = maximum(abs.(imag.(profile[finite_indices])))
        for i in 1:length(profile)
            if !isfinite(profile[i])
                # Maintain the negative imaginary sign for absorption
                profile[i] = -im * max_val
            end
        end
    end
    
    return ManolopoulosCAP(profile)
end
