# API Reference

This page contains the auto-generated documentation for the public functions and types in `SchroedingerEquation.jl`.

```@index
```

## Core Types

```@docs
RealSpaceGrid1D
Wavefunction1D
AbstractBasis
AbstractPotential
AbstractTimeDependentPotential
AbstractHamiltonian
AbstractBoundaryCondition
PeriodicBoundary
HardWallBoundary
AbsorbingBoundary
AbstractWorkspace
Hamiltonian1D
```

## Potentials

```@docs
HarmonicPotential
SquareWellPotential
OscillatingField
CombinedPotential
CustomPotential
CustomTimeDependentPotential
ManolopoulosCAP
```

## Solvers

```@docs
build_hamiltonian
solve_tise
propagate_tdse
propagate_ssfm
step_crank_nicolson
```

## Workspaces

```@docs
TDSEWorkspace
SSFMWorkspace
```

## Utilities

```@docs
strip_units
apply_units
probability_density
expectation_value
build_momentum_operator
LinearAlgebra.normalize!
```

## Visualization

```@docs
animate_evolution
live_evolution
```
