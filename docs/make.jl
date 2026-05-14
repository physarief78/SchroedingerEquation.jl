using Documenter
using SchroedingerEquation
using LinearAlgebra

makedocs(
    sitename = "SchroedingerEquation.jl",
    format = Documenter.HTML(),
    modules = [SchroedingerEquation],
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Basis & Grids" => "basis.md",
            "Hamiltonians" => "hamiltonian.md",
            "Solvers" => "solvers.md",
            "Theory & Methods" => "theory.md",
        ],
        "Python Support" => "python.md",
        "Tutorials" => [
            "Free Space Dispersion" => "tutorials/01_free_space.md",
            "Phase and Observables" => "tutorials/02_observables.md",
            "Quantum Tunneling" => "tutorials/03_tunneling.md",
            "Absorbing Boundaries" => "tutorials/04_absorbing.md",
            "Ultrafast Excitation" => "tutorials/05_excitation.md",
            "Interactive Scattering" => "tutorials/06_interactive.md",
            "Parallel Parameter Sweep" => "tutorials/07_parallel.md",
            "Distributed Ensemble" => "tutorials/08_distributed.md",
        ],
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/physarief78/SchroedingerEquation.jl.git",
    devbranch = "main"
)
