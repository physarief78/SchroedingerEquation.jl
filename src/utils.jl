# Utility functions for wavefunctions

"""
    LinearAlgebra.normalize!(wf::Wavefunction1D)
Normalizes the wavefunction such that integral |psi|^2 dx = 1.
"""
function LinearAlgebra.normalize!(wf::Wavefunction1D)
    dx = wf.basis.dx
    norm_factor = sqrt(sum(abs2, wf.psi) * dx)
    wf.psi ./= norm_factor
    return wf
end

"""
    probability_density(wf::Wavefunction1D)
Returns |psi|^2.
"""
probability_density(wf::Wavefunction1D) = abs2.(wf.psi)

"""
    expectation_value(wf::Wavefunction1D, operator::AbstractMatrix)
Calculates <psi | O | psi>.
"""
function expectation_value(wf::Wavefunction1D, operator::AbstractMatrix)
    dx = wf.basis.dx
    return dot(wf.psi, operator, wf.psi) * dx
end

"""
    expectation_value(wf::Wavefunction1D, f::Function)
Calculates <psi | f(x) | psi>.
"""
function expectation_value(wf::Wavefunction1D, f::Function)
    dx = wf.basis.dx
    x = wf.basis.x
    psi = wf.psi
    return sum(i -> abs2(psi[i]) * f(x[i]), 1:length(x)) * dx
end

# --- Unit Utilities ---

"""
    strip_units(val)
Removes units from a value, converting to atomic units if it's a Quantity.
"""
strip_units(val::Number) = val
strip_units(val::Quantity) = ustrip(auconvert(val))
strip_units(val::AbstractArray) = strip_units.(val)

"""
    apply_units(val, unit)
Adds units back to a value.
"""
apply_units(val, unit) = val * unit
