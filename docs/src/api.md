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
AbstractHamiltonian
AbstractBoundaryCondition
```

## Potentials

```@docs
HarmonicPotential
SquareWellPotential
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
```

## Visualization

```@docs
animate_evolution
live_evolution
```
