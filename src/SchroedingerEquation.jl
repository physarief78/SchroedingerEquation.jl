module SchroedingerEquation

using LinearAlgebra
using SparseArrays
using KrylovKit
using FFTW
using Unitful
using UnitfulAtomic

# --- Abstract Types ---
abstract type AbstractBasis end
abstract type AbstractPotential end
abstract type AbstractTimeDependentPotential <: AbstractPotential end
abstract type AbstractHamiltonian end
abstract type AbstractBoundaryCondition end
abstract type AbstractWorkspace end

# --- Core Structs ---

"""
    RealSpaceGrid1D(x_min, x_max, N)
A uniform 1D grid in real space.
"""
struct RealSpaceGrid1D <: AbstractBasis
    x_min::Float64
    x_max::Float64
    N::Int
    x::Vector{Float64}
    dx::Float64

    function RealSpaceGrid1D(x_min, x_max, N)
        # Handle Unitful inputs by converting to internal Atomic Units
        xm = strip_units(x_min)
        xM = strip_units(x_max)
        x = collect(range(xm, xM, length=N))
        dx = x[2] - x[1]
        new(xm, xM, N, x, dx)
    end
end


"""
    Wavefunction1D(psi, basis)
Represents a 1D wavefunction on a specific basis.
"""
struct Wavefunction1D
    psi::Vector{ComplexF64}
    basis::AbstractBasis
end

# --- Potentials ---

struct CustomPotential <: AbstractPotential
    f::Function
end

(p::CustomPotential)(x) = p.f(x)

struct CustomTimeDependentPotential <: AbstractTimeDependentPotential
    f::Function
end

(p::CustomTimeDependentPotential)(x, t) = p.f(x, t)

# --- Hamiltonian ---

struct Hamiltonian1D{B<:AbstractBasis, P<:AbstractPotential, M<:AbstractMatrix} <: AbstractHamiltonian
    basis::B
    potential::P
    H_matrix::M
end

struct PeriodicBoundary <: AbstractBoundaryCondition end
struct HardWallBoundary <: AbstractBoundaryCondition end

# The Production CAP
struct AbsorbingBoundary <: AbstractBoundaryCondition
    L::Float64        # Length of the absorbing region
    strength::Float64 # Optional strength multiplier (default=1.0 for Manolopoulos optimal)
    
    # The pre-factor based on grid physics
    # We include the strength multiplier for fine-tuning if needed
    # coeff = -im * boundary.strength * (hbar^2 / (2.0 * m * L^2))
end

AbsorbingBoundary(L) = AbsorbingBoundary(L, 1.0)

# --- Workspaces for Parallelism ---

struct TDSEWorkspace{F} <: AbstractWorkspace
    A_fact::F
    B::SparseMatrixCSC{ComplexF64, Int}
    B_psi::Vector{ComplexF64}
end

struct SSFMWorkspace{P, Pi} <: AbstractWorkspace
    P_plan::P
    P_inv_plan::Pi
    U_V::Vector{ComplexF64}
    U_T::Vector{ComplexF64}
end

# --- Includes ---
include("potentials.jl")
include("utils.jl")
include("tise.jl")
include("tdse.jl")
include("visualize.jl") # CairoMakie (Headless-friendly)

# live_evolution is defined here but implemented in the GLMakieExt extension
"""
    live_evolution(...)
Interactive GLMakie visualization. Requires `using GLMakie`.
"""
function live_evolution(args...; kwargs...)
    @warn "live_evolution requires GLMakie. Please run `using GLMakie` first."
end

# --- Exports ---
export AbstractBasis, AbstractPotential, AbstractTimeDependentPotential, AbstractHamiltonian, AbstractBoundaryCondition
export PeriodicBoundary, HardWallBoundary, AbsorbingBoundary
export RealSpaceGrid1D, Wavefunction1D, CustomPotential, CustomTimeDependentPotential, Hamiltonian1D
export HarmonicPotential, SquareWellPotential, OscillatingField, CombinedPotential, ManolopoulosCAP
export probability_density, expectation_value, strip_units, apply_units
export build_hamiltonian, solve_tise, build_momentum_operator
export propagate_tdse, propagate_ssfm, animate_evolution, live_evolution
export TDSEWorkspace, SSFMWorkspace

end # module
