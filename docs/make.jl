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
        ],
        "Python Support" => "python.md",
        "Tutorials" => [
            "1. Free Space Dispersion" => "tutorials/01_free_space.md",
            "2. Phase and Observables" => "tutorials/02_observables.md",
            "3. Quantum Tunneling" => "tutorials/03_tunneling.md",
            "4. Absorbing Boundaries" => "tutorials/04_absorbing.md",
            "5. Ultrafast Excitation" => "tutorials/05_excitation.md",
            "6. Interactive Scattering" => "tutorials/06_interactive.md",
            "7. Parallel Parameter Sweep" => "tutorials/07_parallel.md",
            "8. Distributed Ensemble" => "tutorials/08_distributed.md",
        ],
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/physarief78/SchroedingerEquation.jl.git",
    devbranch = "main"
)
